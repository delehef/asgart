use std::mem;
use std::collections::HashMap;

use structs::ALPHABET;
use divsufsort64::*;
use automaton::Segment;

pub struct Searcher {
    cache: HashMap<u64, (usize, usize)>,
}

static CACHE_LEN: usize = 8;

impl Searcher {
    fn indexize(p: &[u8]) -> u64 {
        // fn value(c: u8) -> u64 {
        // match c {
        //     b'A' => {0}
        //     b'T' => {1}
        //     b'G' => {2}
        //     b'C' => {3}
        //     b'N' => {4}
        //     _    => {panic!("Unknown: {}", c)}
        // }
        // }

        // let mut r: u64 = 0;
        // for i in 0..CACHE_LEN {
        // r = r | (value(p[i]) << (3*i));
        // }
        // r

        unsafe { mem::transmute::<[u8; 8], u64>([p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]]) }
    }

    pub fn new(dna: &[u8], sa: &[idx]) -> Searcher {
        let mut s = Searcher { cache: HashMap::new() };

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
                                                let count = sa_searchb64(dna.as_ptr(),
                                                                         dna.len() as i64,
                                                                         p.as_ptr(),
                                                                         p.len() as i64,
                                                                         sa.as_ptr(),
                                                                         sa.len() as i64,
                                                                         &mut out,
                                                                         0,
                                                                         sa.len() as i64);
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

    pub fn search(&self, dna: &[u8], sa: &[idx], pattern: &[u8]) -> Vec<Segment> {
        let index = Searcher::indexize(pattern);

        let mut out = 0;
        if !self.cache.contains_key(&index) {
            panic!("Unable to find {} ({})",
                   index,
                   String::from_utf8(pattern[0..CACHE_LEN].to_vec()).unwrap());
        }

        let (lstart, rstart) = self.cache[&index];
        let count = unsafe {
            sa_searchb64(dna.as_ptr(),
                         dna.len() as i64,
                         pattern.as_ptr(),
                         pattern.len() as i64,
                         sa.as_ptr(),
                         sa.len() as i64,
                         &mut out,
                         lstart as idx,
                         rstart as idx)
        };

        let mut rr = Vec::new();
        for i in 0..count {
            let start = sa[(out + i) as usize];
            rr.push(Segment {
                tag: 0,
                start: start as usize,
                end: start as usize + pattern.len(),
            });
        }
        rr
    }
}
