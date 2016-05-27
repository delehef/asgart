use std::collections::HashSet;
use std::cmp;
use bio::data_structures::suffix_array::SuffixArray;
use std::fmt;
use bio::alignment::pairwise::*;

const MIN_DUPLICATION_SIZE: usize = 1000;
const MAX_ALIGNMENT_SIZE: usize = 100000;


macro_rules! log(
    ($($arg:tt)*) => (
        // match writeln!(&mut ::std::io::stderr(), $($arg)* ) {
        //     Ok(_) => {},
        //     Err(x) => {} //panic!("Unable to write to stderr: {}", x),
        // }
    )
);


#[derive(Debug, Clone)]
pub struct SD {
    pub left: usize,
    pub right: usize,
    pub size: usize,
    pub rate: f32,
}

enum SearchState {
    Start,
    Grow,
    SparseGrow,
    Proto,
}

#[derive(Debug,Clone)]
pub struct ProtoSD {
    bottom: usize,
    top: usize,
    matches: Vec<Segment>,
}

#[derive(Debug,Clone)]
pub enum ProcessingSD {
    Done(SD),
    ForFuzzy{psd: ProtoSD, strand2_start: usize},
    ForSW{psd: ProtoSD, left_match: Vec<u8>, right_match: Vec<u8>, right_segment: Segment},
    Empty,
}

fn make_right_arms(p: &ProtoSD, max_hole_size: u32) -> Vec<Segment> {
    let matches = p.matches.clone();
    merge_segments_with_delta(&matches, max_hole_size as u64)
}

fn make_duplication(psd: ProtoSD, strand1: &[u8], strand2: &[u8], max_hole_size: u32, align: bool) -> Vec<ProcessingSD> {
    let right_segments = make_right_arms(&psd, max_hole_size);
    let mut r = Vec::new();

    for right_segment in &right_segments {
        if right_segment.end - right_segment.start < MIN_DUPLICATION_SIZE {continue;}

        if align {
            if psd.top - psd.bottom  > MAX_ALIGNMENT_SIZE || right_segment.end - right_segment.start > MAX_ALIGNMENT_SIZE {
                r.push(ProcessingSD::ForFuzzy {
                    psd: psd.clone(),
                    strand2_start: right_segment.start
                });
            } else {
            // Fetch left and right candidate areas
                let left_match = &strand1[psd.bottom..psd.top];
                let right_match = &strand2[right_segment.start..right_segment.end];

                r.push(ProcessingSD::ForSW {
                    psd: psd.clone(),
                    left_match: Vec::from(left_match),
                    right_match: Vec::from(right_match),
                    right_segment: right_segment.clone(),
                })
            }
        } else {
            r.push(ProcessingSD::Done(SD{left: psd.bottom, right: right_segment.start, size: psd.top - psd.bottom, rate: 0.0}));
        }
    }

    r
}

pub fn align_perfect(p: ProcessingSD) -> ProcessingSD {
    if let ProcessingSD::ForSW{psd, left_match, right_match, right_segment} = p {
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(left_match.len(), right_match.len(), -5, -1, &score);
        let alignment = aligner.local(&left_match, &right_match);
        if alignment.operations.len() < MIN_DUPLICATION_SIZE {return ProcessingSD::Empty;};

        let result = SD {
            left: psd.bottom + alignment.xstart,
            right: right_segment.start + alignment.ystart,
            size: alignment.operations.len(),
            rate: 0.0
        };
        println!("Left position:  {}", result.left);
        println!("Right position: {}", result.right);
        println!("Size:           {}", result.size);
        println!("Alignment: \n{}\n", alignment.pretty(&left_match, &right_match));

        ProcessingSD::Done(result)
    } else {
        p
    }
}

pub fn align_fuzzy(strand1: &[u8], strand2: &[u8], p: ProcessingSD) -> ProcessingSD {
    if let ProcessingSD::ForFuzzy{psd, strand2_start} = p {
        let moar = expand_nw(strand1, strand2, psd.bottom, strand2_start);
        ProcessingSD::Done(moar)
    } else {
        p
    }
}

pub fn search_duplications(strand1: &[u8], strand2: &[u8], sa: &SuffixArray, start: usize, end: usize,
                           candidate_size: usize, max_gap_size: u32, align: bool) -> Vec<ProcessingSD> {
    let mut r = Vec::new();

    let mut i = start;
    let mut gap = 0;
    let mut current_start = 0;
    let mut current_segments = Vec::new();
    let mut state = SearchState::Start;

    loop {
        match state {
            SearchState::Start => {
                if i+candidate_size >= end-1 {
                    break;
                }
                gap = 0;
                current_start = i;
                if strand1[i] == b'N' || strand1[i] == b'n' {
                    i += 1;
                    state = SearchState::Start;
                } else {
                    current_segments = Vec::new();
                    let new_segments = search(strand2, sa, &strand1[i..i+candidate_size], candidate_size);
                    if !new_segments.is_empty() {
                        current_segments = merge_segments(&new_segments);
                        state = SearchState::Grow;
                    } else {
                        i += 1;
                        state = SearchState::Start;
                    }
                }

            },
            SearchState::Grow => {
                i += candidate_size/2;
                let mut new_matches = search(strand2, sa, &strand1[i..cmp::min(i+candidate_size, strand1.len()-1)], candidate_size);

                if i >= strand1.len() - candidate_size {
                    state = SearchState::Proto;
                } else if segments_to_segments_distance(&new_matches, &current_segments) <= max_gap_size {
                    current_segments.append(&mut new_matches);
                    current_segments = merge_segments(&current_segments);
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::SparseGrow => {
                i += 1;
                gap += 1;
                let mut new_matches = search(strand2, sa, &strand1[i..i+candidate_size], candidate_size);

                if (gap > max_gap_size) || (i >= strand1.len() - candidate_size) {
                    state = SearchState::Proto;
                } else if segments_to_segments_distance(&new_matches, &current_segments) <= max_gap_size {
                    current_segments.append(&mut new_matches);
                    current_segments = merge_segments(&current_segments);
                    gap = 0;
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::Proto => {
                if i - current_start >= MIN_DUPLICATION_SIZE {
                    let psd = ProtoSD {
                        bottom: current_start,
                        top: i,
                        matches: current_segments.clone()
                    };
                    let mut result = make_duplication(psd, strand1, strand2, max_gap_size, align);
                    r.append(&mut result);
                }
                state = SearchState::Start;
            },
        }
    }

    r
}


pub struct AlignmentResult {
    pub errors: u32,
    pub ins_la: u32,
    pub ins_ra: u32,
}


pub fn needleman_wunsch(q1: &[u8], q2: &[u8]) -> AlignmentResult {
    const DELETE_SCORE: i64 = -5;
    const MISMATCH_SCORE: i64 = -3;
    const MATCH_SCORE: i64 = 1;

    assert!(q1.len() == q2.len());
    let size: usize = q1.len();

    let mut f = vec![0 as i64; size * size];

    // for i in 0..size {
    for (i, fi) in f.iter_mut().enumerate().take(size) {
        *fi = DELETE_SCORE*i as i64;
    }

    for j in 0..size {
        f[j*size as usize] = DELETE_SCORE*j as i64;
    }

    for i in 1..size {
        for j in 1..size {
            let _match = f[(i-1)+(j-1)*size] + if q1[i] == q2[j] {MATCH_SCORE} else {MISMATCH_SCORE};
            let del = f[(i-1)+(j)*size] + DELETE_SCORE;
            let ins = f[i+(j-1)*size] + DELETE_SCORE;
            f[i+j*size] = cmp::max(_match, cmp::max(del, ins));
        }
    }

    let mut i = size-1 as usize;
    let mut j = size-1 as usize;

    let mut result = AlignmentResult {errors: 0, ins_la: 0, ins_ra: 0};
    while i>1 && j>1 {
        if i>0 && j>0
            && (f[i+j*size] == f[(i-1)+(j-1)*size] + if q1[i] == q2[j] {MATCH_SCORE} else {MISMATCH_SCORE})
            {
                if q1[i] != q2[j] {
                    result.errors += 1;
                }
                i -= 1;
                j -= 1;
            } else if i>0 && f[i+j*size] == f[(i-1)+j*size] + DELETE_SCORE {
                result.ins_la += 1;
                i -= 1;
            } else if j>0 && f[i+j*size] == f[i+(j-1)*size] + DELETE_SCORE {
                result.ins_ra += 1;
                j -= 1;
            }
    }
    result
}


pub fn strcmp(s1: &[u8], s2: &[u8]) -> i32 {
    assert!(s1.len() == s2.len());
    let n = s2.len();

    for i in 0..n {
        if s1[i] != s2[i] {
            return s1[i] as i32 - s2[i] as i32;
        }
    }
    0
}


pub fn translate_nucleotide(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _    => panic!("Not a nucleotide: {}", n),
    }
}


pub fn translated(text: &[u8]) -> Vec<u8> {
    text.iter().map(|x| translate_nucleotide(*x)).collect::<Vec<u8>>()
    // let mut r = Vec::with_capacity(text.len());
    // for i in 0..text.len() {
    //     r.push(translate_nucleotide(text[i]));
    // }
    // r
}

#[derive(Clone)]
pub struct Segment {
    start: usize,
    end: usize,
}
impl fmt::Debug for Segment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "S{{{}, {}}} ({})", self.start, self.end, self.end - self.start)
    }
}

fn merge_segments_with_delta(_segments: &[Segment], delta: u64) -> Vec<Segment> {
    let mut segments = _segments.to_vec();
    let mut r = Vec::new();
    segments.sort_by(|x, y| x.start.cmp(&y.start));

    r.push(segments[0].clone());
    for current_segment in segments.iter().skip(1) {
        if current_segment.start - r.last().unwrap().end <= delta as usize{
            r.last_mut().unwrap().end = current_segment.end;
        } else {
            r.push(current_segment.clone());
        }
    }
    r
}

fn merge_segments(_segments: &[Segment]) -> Vec<Segment> {
    let mut segments = _segments.to_vec();
    let mut r = Vec::new();
    segments.sort_by(|x, y| x.start.cmp(&y.start));

    r.push(segments[0].clone());
    for current_segment in segments.iter().skip(1) {
        if r.last().unwrap().end < current_segment.start {
            r.push(current_segment.clone());
        } else if r.last().unwrap().end  < current_segment.end {
            r.last_mut().unwrap().end = current_segment.end;
        }
    }
    r
}

fn segments_to_segments_distance(segments: &[Segment], others: &[Segment]) -> u32 {
    let mut distance = 250000000;
    for other in others {
        for segment in segments {
            if other.start >= segment.start && other.start <= segment.end {
                return 0
            } else {
                let local_distance = cmp::min((segment.start as i64 - other.end as i64).abs(), (segment.end as i64 - other.start as i64).abs());
                if local_distance < distance { distance = local_distance }
            }
        }
    }

    distance as u32
}

pub fn search(dna: &[u8], array: &[usize], pattern: &[u8], candidate_size: usize) -> Vec<Segment> {
    let mut lo = 0;
    let mut hi = dna.len()-1;
    let mut result = HashSet::new();

    while lo <= hi {
        let mid = (lo+hi)/2;
        if array[mid]+pattern.len() > dna.len() {break;}
        let substring = &dna[array[mid]..array[mid]+pattern.len()];
        let res = strcmp(pattern, substring);

        if res > 0 {
            lo = mid + 1;
        } else if res < 0 {
            hi = mid - 1;
        } else { // TIME TO BLOB!
            let mut l = mid;
            let mut r = mid;
            let mut expand_left = true;
            let mut expand_right = true;

            while expand_left || expand_right && array[r]+pattern.len() < dna.len() {
                if expand_left {
                    if *pattern == dna[array[l]..array[l]+pattern.len()-1] {
                        result.insert(array[l]);
                        l -= 1;
                        if l == 0 { expand_left = false; }
                    } else {
                        expand_left = false;
                    }
                }
                if expand_right {
                    if *pattern == dna[array[r]..array[r]+pattern.len()] {
                        result.insert(array[r]);
                        r += 1;
                        if r >= array.len() - 1 { expand_right = false; }
                    } else {
                        expand_right = false;
                    }
                }
            }
            let mut rr = Vec::new();
            for i in result {
                rr.push(Segment{start: i, end: i+candidate_size});
            }
            return rr;
        }
    }
    let mut rr = Vec::new();
    for i in result {
        rr.push(Segment{start: i, end: i+candidate_size});
    }
    rr
}


fn expand_nw(strand1: &[u8], strand2: &[u8], straight_start: usize, reverse_start: usize) -> SD {
    const PRECISION: f32 = 0.8;
    const EXPANSION_STEP: usize = 1000;

    let mut offset:usize = 0;
    let mut rate = 1.0;

    let mut errors = 0;
    let mut correction_la = 0;
    let mut correction_ra = 0;

    let mut la_start;
    let mut ra_start;
    while rate >= PRECISION {
        la_start = straight_start + offset + correction_la;
        ra_start = reverse_start + offset + correction_ra;

        if (la_start+EXPANSION_STEP >= strand1.len()) || (ra_start+EXPANSION_STEP >= strand2.len()) { break; }

        let result = needleman_wunsch(
            &strand1[la_start..la_start + EXPANSION_STEP],
            &strand2[ra_start..ra_start + EXPANSION_STEP]);

        errors += result.errors;

        if result.ins_la > result.ins_ra {
            correction_ra += (result.ins_la - result.ins_ra) as usize;
        } else {
            correction_la += (result.ins_ra - result.ins_la) as usize;
        }

        rate = 1.0 - (errors as f32)/(offset as f32 + EXPANSION_STEP as f32);
        offset += EXPANSION_STEP;
    }

    SD {
        left: straight_start,
        right: reverse_start,
        size: offset+EXPANSION_STEP,
        rate: rate,
    }
}
