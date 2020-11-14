use ::rayon::prelude::*;
use anyhow::{anyhow, Context, Result};
use serde_derive::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use regex::Regex;

pub const COLLAPSED_NAME: &str = "ASGART_COLLAPSED";
pub const ALPHABET: [u8; 5] = [b'A', b'T', b'G', b'C', b'N'];
pub const ALPHABET_MASKED: [u8; 5] = [b'a', b't', b'g', b'c', b'n'];

#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct RunSettings {
    pub probe_size:             usize,
    pub max_gap_size:           u32,
    pub min_duplication_length: usize,
    pub max_cardinality:        usize,
    pub trim:                   Option<(usize, usize)>,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub reverse:     bool,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub complement:  bool,
    pub skip_masked: bool,

    #[serde(skip_serializing)]
    #[serde(default)]
    pub threads_count: usize,
    #[serde(skip_serializing)]
    #[serde(default)]
    pub compute_score: bool,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Start {
    pub name:     String,
    pub position: usize,
    pub length:   usize,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct StrandResult {
    pub name:   String,
    pub length: usize,
    pub map:    Vec<Start>,
}
impl StrandResult {
    pub fn has_chr(&self, name: &str) -> bool {
        self.map.iter().any(|chr| chr.name == name)
    }

    pub fn find_chr(&self, name: &str) -> Option<&Start> {
        self.map.iter().find(|chr| chr.name == name)
    }

    pub fn find_chr_index(&self, name: &str) -> Option<usize> {
        self.map.iter().position(|chr| chr.name == name)
    }

    pub fn find_chr_by_pos(&self, pos: usize) -> Option<&Start> {
        self.map
            .iter()
            .find(|&chr| pos >= chr.position && pos < chr.position + chr.length)
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct RunResult {
    pub strand:   StrandResult,
    pub settings: RunSettings,
    pub families: Vec<SDsFamily>,
}
impl RunResult {
    pub fn from_stdin() -> Result<RunResult> {
        serde_json::from_reader(std::io::stdin())
            .with_context(|| "Failed to parse JSON data from STDIN")
    }

    pub fn from_file(filename: &str) -> Result<RunResult> {
        let mut f = File::open(filename)
            .with_context(|| format!("Cannot read data from `{}`", filename))?;
        let mut s = String::new();
        let _ = f.read_to_string(&mut s);
        serde_json::from_str(&s)
            .with_context(|| format!("Failed to parse JSON data from `{}`", filename))
    }

    pub fn from_files(filenames: &[String]) -> Result<RunResult> {
        let results = filenames
            .iter()
            .map(|filename| RunResult::from_file(filename))
            .collect::<std::result::Result<Vec<RunResult>, _>>()?;

        for result in &results {
            if result.strand.name != results[0].strand.name {
                return Err(anyhow!(
                    "Trying to combine ASGART files from different sources: `{}` and `{}`",
                    result.strand.name,
                    results[0].strand.name,
                ));
            }
        }

        let r = RunResult {
            settings: results[0].settings,
            strand:   results[0].strand.clone(),
            families: results
                .iter()
                .flat_map(|ref r| r.families.iter())
                .cloned()
                .collect::<Vec<_>>(),
        };

        Ok(r)
    }

    pub fn remove_direct(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.reversed));
        self.families.retain(|f| !f.is_empty());
    }

    pub fn remove_reversed(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| !sd.reversed));
        self.families.retain(|f| !f.is_empty());
    }

    pub fn remove_uncomplemented(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.complemented));
        self.families.retain(|f| !f.is_empty());
    }

    pub fn remove_complemented(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| !sd.complemented));
        self.families.retain(|f| !f.is_empty());
    }

    pub fn remove_inter(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.chr_left == sd.chr_right));
        self.families.retain(|f| !f.is_empty());
    }

    pub fn remove_inter_relaxed(&mut self) {
        self.families.iter_mut().for_each(|family| {
            family.retain(|sd| {
                (sd.chr_left == sd.chr_right)
                    || sd.chr_left == COLLAPSED_NAME
                    || sd.chr_right == COLLAPSED_NAME
            })
        });
        self.families.retain(|f| !f.is_empty());
    }

    pub fn remove_intra(&mut self) {
        self.families
            .iter_mut()
            .for_each(|family| family.retain(|sd| sd.chr_left != sd.chr_right));
        self.families.retain(|f| !f.is_empty());
    }

    pub fn max_family_members(&mut self, m: usize) {
        self.families.retain(|family| family.len() <= m)
    }

    pub fn restrict_fragments<T: AsRef<str>>(&mut self, to_keep: &[T]) {
        self.families.iter_mut().for_each(|family| {
            family.retain(|sd| {
                to_keep.iter().any(|n| n.as_ref() == sd.chr_left)
                    && to_keep.iter().any(|n| n.as_ref() == sd.chr_right)
            })
        });
        self.families.retain(|f| !f.is_empty());
        self.strand
            .map
            .retain(|c| to_keep.iter().any(|n| n.as_ref() == c.name));
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position =
                    self.strand.find_chr(&sd.chr_left).unwrap().position + sd.chr_left_position;
                sd.global_right_position =
                    self.strand.find_chr(&sd.chr_right).unwrap().position + sd.chr_right_position;
            }
        }
    }

    pub fn restrict_fragments_regexp<T: AsRef<str>>(&mut self, to_keep: &str) -> anyhow::Result<()> {
        let re = Regex::new(to_keep)?;
        self.families.iter_mut().for_each(|family| {
            family.retain(|sd| {
                re.is_match(&sd.chr_left) && re.is_match(&sd.chr_right)
            })
        });
        self.families.retain(|f| !f.is_empty());
        self.strand
            .map
            .retain(|c| re.is_match(&c.name));
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position =
                    self.strand.find_chr(&sd.chr_left).unwrap().position + sd.chr_left_position;
                sd.global_right_position =
                    self.strand.find_chr(&sd.chr_right).unwrap().position + sd.chr_right_position;
            }
        }

        Ok(())
    }

    pub fn exclude_fragments<T: AsRef<str>>(&mut self, to_exclude: &[T]) {
        self.families.iter_mut().for_each(|family| {
            family.retain(|sd| {
                !to_exclude.iter().any(|n| n.as_ref() == sd.chr_left)
                    && !to_exclude.iter().any(|n| n.as_ref() == sd.chr_right)
            })
        });
        self.families.retain(|f| !f.is_empty());
        self.strand
            .map
            .retain(|c| !to_exclude.iter().any(|n| n.as_ref() == c.name));
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position =
                    self.strand.find_chr(&sd.chr_left).unwrap().position + sd.chr_left_position;
                sd.global_right_position =
                    self.strand.find_chr(&sd.chr_right).unwrap().position + sd.chr_right_position;
            }
        }
    }

    pub fn exclude_fragments_regexp<T: AsRef<str>>(&mut self, to_exclude: &str) -> anyhow::Result<()> {
        let re = Regex::new(to_exclude)?;
        self.families.iter_mut().for_each(|family| {
            family.retain(|sd| {
                !re.is_match(&sd.chr_left) && !re.is_match(&sd.chr_right)
            })
        });
        self.families.retain(|f| !f.is_empty());
        self.strand
            .map
            .retain(|c| !re.is_match(&c.name));
        self.strand.length = self.strand.map.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in self.strand.map.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for f in self.families.iter_mut() {
            for sd in f.iter_mut() {
                sd.global_left_position =
                    self.strand.find_chr(&sd.chr_left).unwrap().position + sd.chr_left_position;
                sd.global_right_position =
                    self.strand.find_chr(&sd.chr_right).unwrap().position + sd.chr_right_position;
            }
        }

        Ok(())
    }

    pub fn flatten(&mut self) {
        if self.strand.map.len() < 2 {
            return;
        }
        let n = self.strand.map.len() as f64;
        let lengths = self
            .strand
            .map
            .iter()
            .map(|c| c.length as f64)
            .collect::<Vec<_>>();
        let avg = lengths.iter().sum::<f64>() / n;
        let std =
            (1.0 / (n - 1.0) * lengths.iter().map(|x| (x - avg).powf(2.0)).sum::<f64>()).sqrt();

        let mut to_flatten = self
            .strand
            .map
            .iter()
            .cloned()
            .filter(|c| c.length as f64 <= avg + std && c.name.len() > 2) // Try not to remove normal but small scaffolds/chromosomes
            .clone()
            .collect::<Vec<_>>();
        let to_flatten_len = to_flatten.iter().map(|c| c.length).sum::<usize>();
        let mut to_keep = self
            .strand
            .map
            .iter()
            .cloned()
            .filter(|c| !to_flatten.iter().any(|r| c.name == r.name))
            .collect::<Vec<_>>();
        let to_keep_len = to_keep.iter().map(|c| c.length).sum::<usize>();

        let mut i = 0;
        for c in to_keep.iter_mut() {
            c.position = i;
            i += c.length;
        }
        for c in to_flatten.iter_mut() {
            c.position = i;
            i += c.length;
        }

        let to_flatten_positions = to_flatten
            .iter()
            .map(|c| (c.name.clone(), c.position))
            .collect::<HashMap<_, _>>();

        self.strand.map = to_keep;
        self.strand.map.push(Start {
            name:     COLLAPSED_NAME.to_string(),
            position: to_keep_len + 1,
            length:   to_flatten_len,
        });

        self.families.par_iter_mut().for_each(|family| {
            family.iter_mut().for_each(|mut sd| {
                let left_match = to_flatten.iter().any(|n| *n.name == sd.chr_left);
                let right_match = to_flatten.iter().any(|n| *n.name == sd.chr_right);
                if left_match {
                    sd.chr_left_position += to_flatten_positions[&sd.chr_left];
                    sd.chr_left = COLLAPSED_NAME.to_string();
                }
                if right_match {
                    sd.chr_right_position += to_flatten_positions[&sd.chr_right];
                    sd.chr_right = COLLAPSED_NAME.to_string();
                }
            })
        });
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProtoSD {
    pub left:  usize,
    pub right: usize,

    pub left_length:  usize,
    pub right_length: usize,

    pub identity:     f32,
    pub reversed:     bool,
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
        let left_arm = &strand[self.left..=self.left + self.left_length];
        let right_arm = &strand[self.right..=self.right + self.right_length];
        let dist = f64::from(bio::alignment::distance::levenshtein(left_arm, right_arm));

        100.0 * (1.0 - dist / (std::cmp::max(self.left_length, self.right_length) as f64))
    }
}
pub type ProtoSDsFamily = Vec<ProtoSD>;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SD {
    pub chr_left:  String,
    pub chr_right: String,

    pub global_left_position:  usize,
    pub global_right_position: usize,

    pub chr_left_position:  usize,
    pub chr_right_position: usize,

    pub left_length:  usize,
    pub right_length: usize,

    #[serde(default)]
    pub left_seq:  Option<String>,
    #[serde(default)]
    pub right_seq: Option<String>,

    pub identity:     f32,
    pub reversed:     bool,
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
