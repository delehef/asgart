use std::collections::HashMap;
use std::mem;

use superslice::Ext;

use crate::automaton::Segment;
use crate::divsufsort::*;
use crate::structs::ALPHABET;

pub struct Searcher {
    cache: HashMap<u64, (usize, usize)>,
    offset: usize,
}

static CACHE_LEN: usize = 8;
const SSE_STRIDE: usize = 16;
const AVX_STRIDE: usize = 32;

/// Like above but with 32 byte slices
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn avx_compare_mask(one: &[u8], two: &[u8]) -> i32 {
    use std::arch::x86_64::*;
    let onev = _mm256_loadu_si256(one.as_ptr() as *const _);
    let twov = _mm256_loadu_si256(two.as_ptr() as *const _);
    let mask = _mm256_cmpeq_epi8(onev, twov);
    !_mm256_movemask_epi8(mask)
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
pub unsafe fn sse_compare_mask(one: &[u8], two: &[u8]) -> i32 {
    use std::arch::x86_64::*;

    // too lazy to figure out the bit-fiddly way to get this mask
    const HIGH_HALF_MASK: u32 = 0b1111_1111_1111_1111_0000_0000_0000_0000;

    debug_assert!(is_x86_feature_detected!("sse4.2"));

    let onev = _mm_loadu_si128(one.as_ptr() as *const _);
    let twov = _mm_loadu_si128(two.as_ptr() as *const _);
    let mask = _mm_cmpeq_epi8(onev, twov);
    (!_mm_movemask_epi8(mask)) ^ HIGH_HALF_MASK as i32
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn ne_idx_avx(one: &[u8], two: &[u8]) -> std::cmp::Ordering {
    let min_len = one.len().min(two.len());
    let mut idx = 0;
    while idx < min_len {
        let stride_len = AVX_STRIDE.min(min_len - idx);
        let mask = avx_compare_mask(
            &one.get_unchecked(idx..idx + stride_len),
            &two.get_unchecked(idx..idx + stride_len),
        );
        // at the end of the slice the mask might include garbage bytes, so
        // we ignore matches that are OOB
        if mask != 0 && idx + (mask.trailing_zeros() as usize) < min_len {
            let i = idx + mask.trailing_zeros() as usize;
            return one[i].cmp(&two[i]);
            // return Some(idx + mask.trailing_zeros() as usize);
        }
        idx += AVX_STRIDE;
    }
    std::cmp::Ordering::Equal
}
#[doc(hidden)]
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
pub unsafe fn ne_idx_sse(one: &[u8], two: &[u8]) -> std::cmp::Ordering {
    let min_len = one.len().min(two.len());
    let mut idx = 0;
    while idx < min_len {
        let stride_len = SSE_STRIDE.min(min_len - idx);
        let mask = sse_compare_mask(
            &one.get_unchecked(idx..idx + stride_len),
            &two.get_unchecked(idx..idx + stride_len),
        );
        if mask != 0 && idx + (mask.trailing_zeros() as usize) < min_len {
            let i = idx + mask.trailing_zeros() as usize;
            return one[i].cmp(&two[i]);
            // return Some(idx + mask.trailing_zeros() as usize);
        }
        idx += SSE_STRIDE;
    }
    std::cmp::Ordering::Equal
}
#[inline]
#[allow(dead_code)]
#[doc(hidden)]
pub fn ne_idx_fallback(one: &[u8], two: &[u8]) -> Option<usize> {
    one.iter().zip(two.iter()).position(|(a, b)| a != b)
}

impl Searcher {
    fn indexize(p: &[u8]) -> u64 {
        unsafe { mem::transmute::<[u8; 8], u64>([p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]) }
    }

    pub fn new(dna: &[u8], sa: &[SAIdx], offset: usize) -> Searcher {
        let mut s = Searcher {
            cache: HashMap::new(),
            offset,
        };

        unsafe {
            for a in &ALPHABET {
                for b in &ALPHABET {
                    for c in &ALPHABET {
                        for d in &ALPHABET {
                            for e in &ALPHABET {
                                for f in &ALPHABET {
                                    for g in &ALPHABET {
                                        for h in &ALPHABET {
                                            let p = vec![*a, *b, *c, *d, *e, *f, *g, *h];
                                            let index = Searcher::indexize(&p);
                                            let (start, count): (usize, usize) = {
                                                let mut out = 0;
                                                let count = sa_searchb64(
                                                    dna.as_ptr(),
                                                    dna.len() as i64,
                                                    p.as_ptr(),
                                                    p.len() as i64,
                                                    sa.as_ptr(),
                                                    sa.len() as i64,
                                                    &mut out,
                                                    0,
                                                    sa.len() as i64,
                                                );
                                                (out as usize, count as usize)
                                            };
                                            s.cache.insert(index, (start, start + count));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        s
    }

    pub fn search(&self, dna: &[u8], sa: &[SAIdx], pattern: &[u8]) -> Vec<Segment> {
        #[inline]
        fn stringcmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
            // return unsafe { ne_idx_sse(a, b) };
            // return unsafe { ne_idx_avx(a, b) };
            a.cmp(b)
        }

        let index = Searcher::indexize(pattern);

        if !self.cache.contains_key(&index) {
            panic!(
                "Unable to find {} ({})",
                index,
                String::from_utf8(pattern[0..CACHE_LEN].to_vec()).unwrap()
            );
        }

        let (lstart, rstart) = self.cache[&index];
        let range = &sa[lstart..rstart].equal_range_by(|x| {
            if *x as usize + pattern.len() > dna.len() {
                std::cmp::Ordering::Less
            } else {
                stringcmp(&dna[*x as usize..*x as usize + pattern.len()], &pattern)
            }
        });

        sa[lstart + range.start..lstart + range.end]
            .iter()
            .map(|start| Segment {
                tag: 0,
                start: self.offset + *start as usize,
                end: self.offset + *start as usize + pattern.len(),
            })
            .collect()
    }
}
