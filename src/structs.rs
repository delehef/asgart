use std::io::prelude::*;
use std::fs::File;
use ::errors::*;
use ::rayon::prelude::*;


pub const COLLAPSED_NAME: &str = "ASGART_COLLAPSED";

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

    pub fn from_files(filenames: &[String]) -> Result<RunResult> {
        let results = filenames
            .iter()
            .map(|filename| RunResult::from_file(filename))
            .collect::<std::result::Result<Vec<RunResult>, _>>()?;

        for result in &results {
            if result.strand.name != results[0].strand.name {
                bail!(format!("Trying to combine ASGART files from different sources: `{}` and `{}`",
                              result.strand.name, results[0].strand.name,
                ));
            }
        }

        let r = RunResult {
            settings: results[0].settings.clone(),
            strand:   results[0].strand.clone(),
            families: results.iter().flat_map(|ref r| r.families.iter()).cloned().collect::<Vec<_>>(),
        };

        Ok(r)
    }

    pub fn remove_direct(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.reversed))
    }

    pub fn remove_reversed(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| !sd.reversed))
    }

    pub fn remove_uncomplemented(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.complemented))
    }

    pub fn remove_complemented(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| !sd.complemented))
    }

    pub fn remove_inter(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.chr_left == sd.chr_right))
    }

    pub fn remove_intra(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.chr_left != sd.chr_right))
    }

    pub fn restrict_fragments<T: AsRef<str>>(&mut self, to_keep: &[T]) {
        self.families.iter_mut().for_each(|family| family.retain(|sd| {
            to_keep.iter().any(|n| n.as_ref() == sd.chr_left)
                && to_keep.iter().any(|n| n.as_ref() == sd.chr_right)
        }));
        self.strand.map.retain(|c| to_keep.iter().any(|n| n.as_ref() == c.name));
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position = self.strand.find_chr(&sd.chr_left).unwrap().position
                    + sd.chr_left_position;
                sd.global_right_position = self.strand.find_chr(&sd.chr_right).unwrap().position
                    + sd.chr_right_position;
            }
        }
    }

    pub fn exclude_fragments<T: AsRef<str>>(&mut self, to_exclude: &[T]) {
        self.families.iter_mut().for_each(|family| family.retain(|sd| {
            !(to_exclude.iter().any(|n| n.as_ref() == sd.chr_left)
              || to_exclude.iter().any(|n| n.as_ref() == sd.chr_right))
        }));
        self.strand.map.retain(|c| !to_exclude.iter().any(|n| n.as_ref() == c.name));
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position = self.strand.find_chr(&sd.chr_left).unwrap().position
                    + sd.chr_left_position;
                sd.global_right_position = self.strand.find_chr(&sd.chr_right).unwrap().position
                    + sd.chr_right_position;
            }
        }
    }

    pub fn flatten(&mut self) {
        if self.strand.map.len() < 2 { return; }
        let n = self.strand.map.len() as f64;
        let lengths = self.strand.map.iter().map(|c| c.length as f64).collect::<Vec<_>>();
        let avg = lengths.iter().sum::<f64>()/n;
        let std = (1.0/(n - 1.0) * lengths.iter().map(|x| (x - avg).powf(2.0)).sum::<f64>()).sqrt();

        let to_remove = self.strand.map
            .iter().cloned()
            .filter(|c| c.length as f64 <= avg + std)
            .map(|c| c.name.clone())
            .collect::<Vec<String>>();

        self.families
            .par_iter_mut()
            .for_each(|family|
                      family.iter_mut()
                      .for_each(|mut sd| {
                          let left_match = to_remove.iter().any(|n| *n == sd.chr_left);
                          let right_match = to_remove.iter().any(|n| *n == sd.chr_right);
                          if left_match {
                              sd.chr_left = COLLAPSED_NAME.to_string();
                              sd.global_left_position = 0;
                              sd.chr_left_position = 0;
                              sd.left_length = 0;
                          }
                          if right_match {
                              sd.chr_right = COLLAPSED_NAME.to_string();
                              sd.global_right_position = 0;
                              sd.chr_right_position = 0;
                              sd.right_length = 0;
                          }}));

        self.strand.map.retain(|c| !to_remove.iter().any(|n| *n == c.name));
        self.strand.map.sort_by(|a, b| a.name.cmp(&b.name));
        self.strand.map.push(Start{name: COLLAPSED_NAME.to_string(), position: 0, length: 0});
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();
        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position = self.strand.find_chr(&sd.chr_left).unwrap().position
                    + sd.chr_left_position;
                sd.global_right_position = self.strand.find_chr(&sd.chr_right).unwrap().position
                    + sd.chr_right_position;
            }
        }
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
