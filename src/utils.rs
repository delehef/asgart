use std::collections::HashSet;
use std::cmp;
use std::fmt;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;

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

impl SD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.left, self.size)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.right, self.size)
    }
}

#[derive(Clone)]
pub struct Segment {
    tag: usize,
    start: usize,
    end: usize,
}

impl fmt::Debug for Segment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "S{{{}, {}}} ({})", self.start, self.end, self.end - self.start)
    }
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

fn make_duplications(psd: ProtoSD, strand1: &[u8], strand2: &[u8], max_hole_size: u32, align: bool) -> Vec<ProcessingSD> {
    let right_segments = make_right_arms(&psd, max_hole_size);
    let mut r = Vec::new();

    for right_segment in right_segments {

        // Ignore too small SD
        // let size = cmp::min(psd.top - psd.bottom, right_segment.end - right_segment.start);
        let size = right_segment.end - right_segment.start;
        if size < MIN_DUPLICATION_SIZE {continue;}

        // Ignore N-dominant SD
        if *(&strand1[psd.bottom..psd.top].iter().filter(|c| **c == b'n' || **c == b'N').count()) as f32> 0.1*(psd.top - psd.bottom) as f32 {continue;}

        if ((right_segment.start as i32 - psd.bottom as i32).abs() <= (psd.top - psd.bottom) as i32)
            || (right_segment.start == psd.bottom)
            {continue;}

        if align {
            let ppsd = if right_segment.start != 0 {
                ProtoSD {
                    bottom: right_segment.tag,
                    top: psd.top,
                    matches: psd.matches.clone(),
                }
            } else {
                psd.clone()
            };

            if size > MAX_ALIGNMENT_SIZE {
                r.push(ProcessingSD::ForFuzzy {
                    psd: ppsd,
                    strand2_start: right_segment.start
                });
            } else {
                let left_match = &strand1[psd.bottom..psd.bottom+size];
                let right_match = &strand2[right_segment.start..right_segment.start+size];

                r.push(ProcessingSD::ForSW {
                    psd: psd.clone(),
                    left_match: Vec::from(left_match),
                    right_match: Vec::from(right_match),
                    right_segment: right_segment.clone(),
                })
            }
        } else {
            r.push(ProcessingSD::Done(SD{
                left: cmp::max(psd.bottom, right_segment.tag),
                right: right_segment.start,
                size: size,
                rate: 0.0
            }));
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

pub fn search_duplications(strand1: &[u8], strand2: &[u8], sa: &[usize], start: usize, end: usize,
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
                    current_segments = search(strand2, sa, &strand1[i..i+candidate_size], candidate_size)
                        .into_iter().filter(|x| x.start != i).collect();
                    if current_segments.is_empty() {
                        i += 1;
                        state = SearchState::Start;
                    } else {
                        state = SearchState::Grow;
                    }
                }

            },
            SearchState::Grow => {
                i += candidate_size;
                if strand1[i] == b'N' || strand1[i] == b'n' {
                    state = SearchState::SparseGrow;
                } else {
                    let new_matches: Vec<Segment> = search(strand2, sa, &strand1[i..cmp::min(i+candidate_size, strand1.len()-1)], candidate_size)
                        .into_iter().filter(|x| x.start != i).collect();

                    if i >= strand1.len() - candidate_size {
                        state = SearchState::Proto;
                    } else if segments_to_segments_distance(&new_matches, &current_segments) <= max_gap_size {
                        // current_segments.append(&mut new_matches);
                        // current_segments = merge_segments(&current_segments);
                        append_merge_segments(&mut current_segments, &new_matches, max_gap_size, i);
                        // merge_or_drop_segments(&mut current_segments, &new_matches, max_gap_size);

                        state = SearchState::Grow;
                    } else {
                        state = SearchState::SparseGrow;
                    }
                }
            },
            SearchState::SparseGrow => {
                i += 1;
                gap += 1;
                if (gap > max_gap_size) || (i >= strand1.len() - candidate_size) {
                    state = SearchState::Proto;
                } else if strand1[i] != b'N' && strand1[i] != b'n' {
                    let new_matches = search(strand2, sa, &strand1[i..i+candidate_size], candidate_size);
                    if segments_to_segments_distance(&new_matches, &current_segments) <= max_gap_size {
                        // current_segments.append(&mut new_matches);
                        // current_segments = merge_segments(&current_segments);
                        append_merge_segments(&mut current_segments, &new_matches, max_gap_size, i);
                        // merge_or_drop_segments(&mut current_segments, &new_matches, max_gap_size);
                        gap = 0;
                        state = SearchState::Grow;
                    } else {
                        state = SearchState::SparseGrow;
                    }
                }
            },
            SearchState::Proto => {
                if i - current_start >= MIN_DUPLICATION_SIZE {
                    let psd = ProtoSD {
                        bottom: current_start,
                        top: i,
                        matches: current_segments.clone()
                    };
                    let mut result = make_duplications(psd, strand1, strand2, max_gap_size, align);
                    r.append(&mut result);
                }
                state = SearchState::Start;
            },
        }
    }

    r
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
        b'N' | _ => b'N',
    }
}

pub fn translated(text: &[u8]) -> Vec<u8> {
    text.iter().map(|x| translate_nucleotide(*x)).collect::<Vec<u8>>()
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

fn append_merge_segments(_originals: &mut Vec<Segment>, _news: &[Segment], delta: u32, i: usize) {
    for n in _news {
        let mut added = false;
        for o in _originals.iter_mut() {
            if n.start >= o.start && n.start <= (o.end + delta as usize) &&
                n.end > (o.end + delta as usize) {
                o.end = n.end;
                added = true;
                continue;
            }
        }
        if !added {
            _originals.push(Segment{tag: i, start: n.start, end: n.end});
        }
    }
}

fn merge_or_drop_segments(_originals: &mut Vec<Segment>, _news: &[Segment], delta: u32) {
    for n in _news {
        for o in _originals.iter_mut() {
            if n.start >= o.start && n.start <= (o.end + delta as usize) &&
                n.end > (o.end + delta as usize) {
                o.end = n.end;
                continue;
            }
        }
    }
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
                rr.push(Segment{tag: 0, start: i, end: i+candidate_size});
            }
            return rr;
        }
    }
    let mut rr = Vec::new();
    for i in result {
        rr.push(Segment{tag: 0, start: i, end: i+candidate_size});
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

        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let left = &strand1[la_start..la_start + EXPANSION_STEP];
        let right = &strand2[ra_start..ra_start + EXPANSION_STEP];
        let mut aligner = Aligner::with_capacity(left.len(), right.len(), -5, -1, &score);
        let alignment = aligner.local(left, right);
        errors += alignment.operations.iter().filter(|x| **x == AlignmentOperation::Subst).count();

        correction_la += alignment.operations.iter().filter(|x| **x == AlignmentOperation::Del).count();
        correction_ra += alignment.operations.iter().filter(|x| **x == AlignmentOperation::Ins).count();

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
