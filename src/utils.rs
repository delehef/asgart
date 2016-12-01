pub fn translate_nucleotide(n: u8) -> u8 {
    match n {
        b'A'     => b'T',
        b'T'     => b'A',
        b'G'     => b'C',
        b'C'     => b'G',
        b'N' | _ => b'N',
    }
}

pub fn translated(text: &[u8]) -> Vec<u8> {
    text.iter().map(|x| translate_nucleotide(*x)).collect::<Vec<u8>>()
}
