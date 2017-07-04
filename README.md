# ASGART: a large duplications finder

Asgart (A Segmental duplications Gathering and Refinement Tool) is a tool
designed to search for large duplications amongst one or two DNA strands.

## Licensing

Asgart is distributed under the GPLv3 license. Please see the LICENSE file.

# Main features

# Installation

## Linux

Static binaries for Linux are available for [x86]() and [x86-64]() platforms.

## MacOS

Binaries for macOS are available [here]().

## Windows

Binaries for Windows are available [here]().

## From sources

To build ASGART from sources, you need CMake, a C compiler and the Rust compiler.

Once these requirements are satisfied, you can build ASGART with

```
cargo build --release
```

Once the build finished, you'll find the binary in `target/release/asgart`.


# Usage

## Simple usage

```
asgart seq.fasta seq.fasta 20 100
```

This command will look for duplications in the `seq.fasta` file then write them in a
JSON file in the folder where it was launched.

## Options

# Plotting

# Formats

# Documentation
