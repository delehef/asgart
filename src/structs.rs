#[derive(RustcEncodable, RustcDecodable)]
pub struct Strand {
    pub name: String,
    pub length: usize,
    pub reversed: bool,
    pub translated: bool,
}


#[derive(RustcEncodable, RustcDecodable)]
pub struct RunResult {
    pub strand1: Strand,
    pub strand2: Strand,

    pub sds: Vec<SD>,
}


#[derive(Debug, Clone, RustcEncodable, RustcDecodable)]
pub struct SD {
    pub left: usize,
    pub right: usize,
    pub size: usize,
    pub identity: f32,
}

impl SD {
    pub fn left_part(&self) -> (usize, usize) {
        (self.left, self.size)
    }

    pub fn right_part(&self) -> (usize, usize) {
        (self.right, self.size)
    }
}
