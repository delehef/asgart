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

        _    => b'N'
    }
}

pub fn complemented(text: &[u8]) -> Vec<u8> {
    text.iter().map(|x| complement_nucleotide(*x)).collect::<Vec<u8>>()
}
