pub fn complement_nucleotide(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',

        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        b'n' => b'n',

        _ => b'N',
    }
}

pub fn complemented(text: &[u8]) -> Vec<u8> {
    text.iter()
        .map(|x| complement_nucleotide(*x))
        .collect::<Vec<u8>>()
}

pub fn slugify(x: &str) -> String {
    x.trim()
        .replace([' ', ':', '|'], "_")
}

pub fn make_out_filename(
    filename: Option<&str>,
    default: &str,
    extension: &str,
) -> std::path::PathBuf {
    filename
        .map(|f| {
            let mut path = std::path::PathBuf::from(f);
            if path.is_dir() {
                path.push(default);
            }
            path.set_extension(extension);
            path
        })
        .unwrap_or_else(|| {
            let mut new = std::path::PathBuf::from(default);
            new.set_extension(extension);
            new
        })
}
