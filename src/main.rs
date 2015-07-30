use std::collections::HashMap;

#[derive(Debug)]
struct Node {
    children: HashMap<char, Box<Node>>,
    indexes: Vec<u64>
}


impl Node {
    fn new() -> Node {
        Node {
            children: HashMap::new(),
            indexes: Vec::new()
        }
    }


    fn insert_suffix(&mut self, s: &str, index: u64) {
        self.indexes.push(index);
        if s.len() > 0 {
            let c_index = s.chars().nth(0).unwrap();

            if !self.children.contains_key(&c_index) {
                self.children.insert(c_index, Box::new(Node::new()));
            }

            self.children.get_mut(&c_index).unwrap().insert_suffix(&s[1..], index+1);
        }
    }


    fn search(&self, s: &str) -> Option<Vec<u64>> {
        if s.len() == 0 {
            println!("{:?}", self.children);
            return Some(self.indexes.clone())
        }

        let next_char = &s.chars().nth(0).unwrap();
        if self.children.contains_key(next_char) {
            return self.children[next_char].search(&s[1..])
        }

        None
    }
}

struct SuffixTrie {
    root: Node
}


impl SuffixTrie {
    pub fn new(text: &str) -> SuffixTrie {
        let mut trie = SuffixTrie {root: Node::new()};
        for i in 0..text.chars().count() {
            trie.root.insert_suffix(&text[i..], i as u64);
            println!("{} {}", i, &text[i..]);
        }
        trie
    }

    pub fn search(&self, s: &str) {
        let mut result = self.root.search(s);
        match result {
            None => {println!("Pattern not found")}
            Some(r) => {
                for i in r {
                    println!("{} found at {}", s, i)
                }
            }
        }
    }
}



fn main() {
    let test = "geeksforgeeks.org";
    let tree = SuffixTrie::new(test);

    for i in vec!["ee", "geek", "quiz", "forgeeks"] {
        println!("Looking for {}", i);
        tree.search(i);
    }
}
