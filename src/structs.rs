pub const ALPHABET: [u8; 5] = [b'A', b'T', b'G', b'C', b'N'];

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
}


#[derive(RustcEncodable, RustcDecodable)]
pub struct RunResult {
    pub strand1: StrandResult,
    pub strand2: StrandResult,

    pub kmer: usize,
    pub gap: usize,

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
