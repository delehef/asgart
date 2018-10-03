pub const ALPHABET: [u8; 5] = [
    b'A', b'T', b'G', b'C', b'N'
];
pub const ALPHABET_MASKED: [u8; 5] = [
    b'a', b't', b'g', b'c', b'n'
];

#[derive(Serialize, Deserialize, Clone, Copy)]
pub struct RunSettings {
    pub probe_size: usize,
    pub max_gap_size: u32,
    pub min_duplication_length: usize,
    pub max_cardinality: usize,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub reverse: bool,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub complement: bool,
    pub skip_masked: bool,

    pub interlaced: bool,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub start: usize,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub end: usize,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Start {
    pub name: String,
    pub position: usize,
    pub length: usize,
}

#[derive(Serialize, Deserialize)]
pub struct StrandResult {
    pub name: String,
    pub length: usize,
    pub map: Vec<Start>,
}

impl StrandResult {
    pub fn has_chr(&self, name: &str) -> bool {
        self.map.iter().any(|chr| chr.name == name)
    }

    pub fn find_chr(&self, name: &str) -> &Start {
        self.map.iter().find(|chr| chr.name == name).unwrap()
    }

    pub fn find_chr_by_pos(&self, pos: usize) -> &Start {
        self.map.iter().find(|&chr| pos> chr.position &&  pos < chr.position + chr.length).unwrap()
    }

    pub fn find_chr_index(&self, pos: usize) -> Option<usize> {
        self.map.iter().position(|ref chr| pos> chr.position &&  pos < chr.position + chr.length)
    }
}


#[derive(Serialize, Deserialize)]
pub struct RunResult {
    pub strand1: StrandResult,
    pub strand2: StrandResult,
    pub settings: RunSettings,
    pub sds: Vec<SD>,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SD {
    pub left: usize,
    pub right: usize,
    pub length: usize,
    pub identity: f32,
    pub reversed: bool,
    pub complemented: bool,
}

impl SD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.left, self.length)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.right, self.length)
    }
}
