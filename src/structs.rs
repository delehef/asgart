pub const ALPHABET: [u8; 5] = [
    b'A', b'T', b'G', b'C', b'N'
];
pub const ALPHABET_MASKED: [u8; 5] = [
    b'a', b't', b'g', b'c', b'n'
];

#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct RunSettings {
    pub probe_size:             usize,
    pub max_gap_size:           u32,
    pub min_duplication_length: usize,
    pub max_cardinality:        usize,
    pub trim:                   Option<(usize, usize)>,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub reverse:                bool,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub complement:             bool,
    pub skip_masked:            bool,

    pub interlaced:             bool,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub start:                  usize,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub end:                    usize,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub threads_count:          usize,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub chunk_size:             usize,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub compute_score:          bool,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Start {
    pub name: String,
    pub position: usize,
    pub length: usize,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct StrandResult {
    pub name: String,
    pub length: usize,
    pub map: Vec<Start>,
}

impl StrandResult {
    pub fn has_chr(&self, name: &str) -> bool {
        self.map.iter().any(|chr| chr.name == name)
    }

    // TODO Actually use an error type like civilised humans
    pub fn find_chr(&self, name: &str) -> &Start {
        self.map.iter().find(|chr| chr.name == name).unwrap()
    }

    pub fn find_chr_index(&self, name: &str) -> Option<usize> {
        self.map.iter().position(|chr| chr.name == name)
    }

    pub fn find_chr_by_pos(&self, pos: usize) -> &Start {
        self.map.iter().find(|&chr| pos> chr.position &&  pos < chr.position + chr.length).unwrap()
    }
}


#[derive(Serialize, Deserialize, Debug)]
pub struct RunResult {
    pub strand1: StrandResult,
    pub strand2: StrandResult,
    pub settings: RunSettings,
    pub sds: Vec<SD>,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProtoSD {
    pub left: usize,
    pub right: usize,
    pub length: usize,
    pub identity: f32,
    pub reversed: bool,
    pub complemented: bool,
}
impl ProtoSD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.left, self.length)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.right, self.length)
    }

    pub fn levenshtein(&self, trim: usize, strand1: &[u8], strand2: &[u8]) -> f64 {
        let left_arm  = &strand1[self.left  ..= self.left + self.length];
        let right_arm = &strand2[self.right - trim ..= self.right - trim + self.length];
        let dist = bio::alignment::distance::levenshtein(left_arm, right_arm) as f64;

        let score = 100.0 * (1.0 - dist/(self.length as f64));

        score
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SD {
    pub chr_left: String,
    pub chr_right: String,

    pub global_left_position: usize,
    pub global_right_position: usize,

    pub chr_left_position: usize,
    pub chr_right_position: usize,

    pub length: usize,
    pub identity: f32,
    pub reversed: bool,
    pub complemented: bool,
}
impl SD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.global_left_position, self.length)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.global_right_position, self.length)
    }
}
