pub const ALPHABET: [u8; 5] = [b'A', b'T', b'G', b'C', b'N'];

#[derive(RustcEncodable, RustcDecodable, Clone, Copy)]
pub struct RunSettings {
    pub probe_size: usize,
    pub max_gap_size: u32,
    pub min_duplication_length: usize,
    pub max_cardinality: usize,

    pub reverse: bool,
    pub translate: bool,
    pub interlaced: bool,
    pub skip_masked: bool,
    pub start: usize,
    pub end: usize,
}

#[derive(RustcEncodable, RustcDecodable, Clone)]
pub struct Start {
    pub name: String,
    pub position: usize,
    pub length: usize,
}

#[derive(RustcEncodable, RustcDecodable)]
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


#[derive(RustcEncodable, RustcDecodable)]
pub struct RunResult {
    pub strand1: StrandResult,
    pub strand2: StrandResult,
    pub settings: RunSettings,
    pub sds: Vec<SD>,
}


#[derive(Debug, Clone, RustcEncodable, RustcDecodable)]
pub struct SD {
    pub left: usize,
    pub right: usize,
    pub length: usize,
    pub identity: f32,
    pub reversed: bool,
    pub translated: bool,
}

impl SD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.left, self.length)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.right, self.length)
    }
}
