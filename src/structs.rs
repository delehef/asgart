use std::io::prelude::*;
use std::fs::File;
use ::errors::*;



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

    #[serde(skip_serializing)]
    #[serde(default)]
    pub threads_count:          usize,
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

    pub fn find_chr(&self, name: &str) -> Option<&Start>{
        self.map.iter().find(|chr| chr.name == name)
    }

    pub fn find_chr_index(&self, name: &str) -> Option<usize> {
        self.map.iter().position(|chr| chr.name == name)
    }

    pub fn find_chr_by_pos(&self, pos: usize) -> Option<&Start> {
        self.map.iter().find(|&chr| pos >= chr.position &&  pos < chr.position + chr.length)
    }
}


#[derive(Serialize, Deserialize, Debug)]
pub struct RunResult {
    pub strand: StrandResult,
    pub settings: RunSettings,
    pub families: Vec<SDsFamily>,
}
impl RunResult {
    pub fn from_file(filename: &str) -> Result<RunResult> {
        let mut f = File::open(filename).chain_err(|| format!("Unable to open {}", filename))?;
        let mut s = String::new();
        let _ = f.read_to_string(&mut s);
        serde_json::from_str(&s).chain_err(|| "Failed to parse JSON")
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProtoSD {
    pub left: usize,
    pub right: usize,

    pub left_length: usize,
    pub right_length: usize,

    pub identity: f32,
    pub reversed: bool,
    pub complemented: bool,
}
impl ProtoSD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.left, self.left_length)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.right, self.right_length)
    }

    pub fn levenshtein(&self, strand: &[u8]) -> f64 {
        // TODO take into account R/C duplications
        let left_arm  = &strand[self.left  ..= self.left + self.left_length];
        let right_arm = &strand[self.right ..= self.right + self.right_length];
        let dist = f64::from(bio::alignment::distance::levenshtein(left_arm, right_arm));

        100.0 * (1.0 - dist/(std::cmp::max(self.left_length, self.right_length) as f64))
    }
}
pub type ProtoSDsFamily = Vec<ProtoSD>;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SD {
    pub chr_left: String,
    pub chr_right: String,

    pub global_left_position: usize,
    pub global_right_position: usize,

    pub chr_left_position: usize,
    pub chr_right_position: usize,

    pub left_length: usize,
    pub right_length: usize,

    pub identity: f32,
    pub reversed: bool,
    pub complemented: bool,
}
impl SD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.global_left_position, self.left_length)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.global_right_position, self.right_length)
    }
}
pub type SDsFamily = Vec<SD>;
