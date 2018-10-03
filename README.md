# ASGART: a large duplications finder

`asgart` (A Segmental duplications Gathering and Refinement Tool) is a
multiplatform (GNU/Linux, macOS, Windows) tool designed to search for
large duplications amongst one or two DNA strands.


## Licensing

Asgart is distributed under the GPLv3 license. Please see the LICENSE
file.

# Why should I use ASGART?

![A map of the Human genome long segmental
duplications](screenshots/chord.png)

You should use ASGART if

- you want to find segmental duplications, either direct, reversed
  and/or complement in a DNA sequence;

- you want to find highly similar parts inbetween sequences up to the
  genome scale;

- you want to map highly similar sequences amongst genomes;

- you need an easy way to visualize the results.

# Installation

## Linux

Static binaries for Linux are available [here](https://github.com/delehef/asgart/releases) for x86_64 platforms.

## MacOS

Binaries for macOS are available [here](https://github.com/delehef/asgart/releases).

## Windows

Binaries for Windows are not yet available.

## From sources

To build ASGART from sources, you need CMake, a C compiler and the
[Rust compiler](https://www.rust-lang.org/en-US/install.html).

Once these requirement are installed, clone the repository

```
git clone https://github.com/delehef/asgart.git
cd asgart
git submodule init
git submodule update
```

You can then build ASGART by running the Rust building tool

```
cargo build --release
```

Once the build is finished, you will find the binary in `target/release/`.


# Usage

## Simple usage

First, let us take a look at a simple example:

```
asgart seq.fasta seq.fasta
```

This command will look for duplications in the `seq.fasta` file, then
write them in a JSON file in the folder where it was launched. ASGART
will probe using 20-mers, and guarantee that no duplication will
include gaps longer than 100bp in their arm-to-arm pairwise alignment.

If you wish to look for reversed-complemented duplications, use the
`-R` and `-C` options, that can be combined in `-RC`. And the `-v`
option will give you more informations, as well as a visual overview
of the progress.

```
asgart seq.fasta seq.fasta -RCv
```

## Input

As input, ASGART takes FASTA files containing the sequences within
which to look for duplications. They can be either in the FASTA or
multiFASTA format.

## Output

### JSON

By default, ASGART will write its result in a JSON file in the folder
where it was launched, following the following structure:

```
{
        "strand1": {
                "name": first strand filename,
                "length": FATA file length,
                "map": [
                        {
                                "name": FASTA fragment name,
                                "position": offset in the FASTA file,
                                "length": FASTA fragment length
                        }
                ]
        },

        "strand2": {
                "name": second strand filename,
                "length": FATA file length,
                "map": [
                        {
                                "name": FASTA fragment name,
                                "position": offset in the FASTA file,
                                "length": FASTA fragment length
                        }
                ]
        },

        "settings": {
                "probe_size": probe size used,
                "max_gap_size": maximal gap size used,
                "min_duplication_length": minimal length for a duplicon,
                "max_cardinality": maximal size of a family,
                "skip_masked": were masked nucleotides skipped?,
                "interlaced": were interlaced looked for?
        },

        "sds": [
                {
                        "left": position of the left arm in the first file,
                        "right": position of the right arm in the second file,
                        "length": length of the duplication (bp),
                        "reversed": true if the duplication is reversed, false else,
                        "complemented": true if the duplication is complemented, false else
                },
                ...
        ]
}
```

### GFF

ASGART can also write its results in GFF2 or GFF3 files by using the
`--format` option. For instance, use `--format gff3` to save the
results in a GFF3 file.

## Options

### Functional

  - `--verbose`/`-v` display mnore information and a progress bar

  - `--reverse`/`-R` look for duplication which second arm is reversed

  - `--complement`/`-C` look for duplication which second arm is
    complemented

  - `--max-cardinality` specifies the maximal count of members in a
    duplication family (default: 1000)

  - `--min-length SIZE` specifies the minimal length (in bp) over
    which a duplication is kept in the final result and not discarded
    (default: 1000)

  - `--skip-masked`/`-S` skip soft-masked zones, _i.e._ lowercased
    parts of the input files (default: no)

### Technical

  - `-h`, `--help` display an help screen

  - `--out FILENAME` specifies the file in which the results will be
    written

  - `--prefix NAME` defines a prefix to prepend to the standard out
    file name

  - `--format OUT_FORMAT` sets the output format. Default is `json`,
    but can be set to gff2 or gff3

  - `--threads COUNT` set the numbers of thread to use. Defaults to
    the number of cores abailable on the CPU

  - `--trim START END` run ASGART only on the specified area of the
    first file

# Plotting

ASGART comes with a plotting tool, producing a visual overview of the
duplications. Currently, two type of graphs are available: chord
graphs, or flat graphs.

## Options

  - `-h`, `--help` display an help screen

  - `--out FILENAME` set output file name

  - `--min-length` set the minimal length (in bp) for a duplication to
    be plotted (default: 5000bp)

  - `--min-identity` set the minimal identity rate (in %) for a
    duplication to be plotted (default: 0%).

  - `--no-direct` do not plot direct duplications

  - `--no-reversed` do not plot reversed duplications

  - `--no-uncomplemented` do not plot non-complemented duplications

  - `--no-complemented` do not plot complemented duplications

  - `--features FILE` add an additional track containing features to
    plot alongside the duplications.

  - `--filter-features DISTANCE` don't plot duplications that are
    farther away then `DISTANCE` bp from the features in the track.

### Feature file format

The feature file format contains a list of lines with three values
separated by semi-colons.

1. The label of the feature.
2. the start of the feaure. It may either be a single integer
   representing its absolute coordinate, or be of the form
   `NAME+OFFSET`, defining a start position at `OFFSET` from the start
   of `NAME` chromosomes (from the input FASTA file).
3. The length of the feaure in base pairs.

Comment lines starts with a `#`.

#### Example

```
# This is a comment line
# This is a feature named MYH14, 122358bp long, and starting at the 50,188,186th base of the chromosome 19
MYH14;19+50188186;122358
# This is a feature named Foo, starting on the 123,456,789th base of the input FASTA file and 1250bp long
Foo;123456789;1250
```

## Chord graphs

A chord graph represents duplications amongst a DNA fragment as arcs
linking point on a circle figuring a fragment bend over itself. Their
width is directly proportional to the length of the duplications they
represent.

### Example

`asgart-plot human_genome.json chord --out=flat.svg --min-length 20000`

![Chord graph example](screenshots/chord.png)

## Flat graphs

Flat graphs are made of two superposed horizontal lines, representing
the two fragments analyzed by ASGART, with lines linking left and
right parts of the duplications found, their width proportional to the
length of the duplication.

### Example

`asgart-plot human_Y.json flat --out=flat.svg --no-direct --no-uncomplemented --min-length 2000`

![Flat graph example](screenshots/flat.png)

# Update log

## v1.3.1

- Fix erreneous GFF3 output: seq names are now corrent, no superfluous underscore and correct, relative positions instead of absolute ones.

## v1.3

- Add a new plot format, _genomic_
- Relabel “translate” to “complement”
- Fix the lack of color in SVG export
- `asgart-plot` can now read features tracks, either in custom or GFF3 format
- Add a setting to skip soft-masked zones
- Update dependencies

## v1.2

- Deep refactoring of the plotting system

## v1.1

- Add GFF2 & GFF3 export formats
- Improve build system
- Refactoring
- Fix various small bugs

## v1.0

- First published version
