# Asgart

Asgart (A Segmental duplications Gathering and Refinement Tool) is a tool
designed to computationally search for segmental duplications in a DNA string.

Its key feature is its adaptability to the available computing power by its ability
to comprommise between the quality and precision of the results and the compute time
and memory required.

Asgart is distributed under the GPLv3 license.

# Build and run

You may run Asgart whether in `debug` or `release` with the following command:

```
cargo run [--release] KMER_SIZE MAX_HOLE_SIZE
```
