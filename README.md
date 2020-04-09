# ASGART: A Large Duplications Finder

[ASGART is a multiplatform (GNU/Linux, macOS, Windows)](https://academic.oup.com/bioinformatics/article/34/16/2708/4948616), efficient, tool designed to search for large (>1000bp) duplications amongst one or more DNA sequences, up to the genome scale.


## Licensing

Asgart is distributed under the GPLv3 license. Please see the LICENSE
file.

[How to cite?](https://scholar.google.com/scholar?hl=fr&as_sdt=0%2C5&q=asgart&btnG=#d=gs_cit&u=%2Fscholar%3Fq%3Dinfo%3AHI7UCO1nHU8J%3Ascholar.google.com%2F%26output%3Dcite%26scirp%3D2%26hl%3Dfr)

# Why Should I Use ASGART?

![A map of the Human genome long segmental duplications](screenshots/chord.png)

You should use ASGART if

- you want to find segmental duplications in one or more DNA sequences;

- you want to find highly similar parts between sequences up to the
  genome scale;

- you want to map highly similar sequences amongst genomes;

- you need an easy way to visualize the results.

# Installation

## Linux

Static binaries for Linux are available [here](https://github.com/delehef/asgart/releases) for x86_64 platforms.

## MacOS

Binaries for macOS are available [here](https://github.com/delehef/asgart/releases).

## From Sources

To build ASGART from sources, you need CMake, a C compiler and the
[Rust compiler](https://www.rust-lang.org/en-US/install.html).

Once these requirement are satisfied, clone the repository and its submodule

```
git clone https://github.com/delehef/asgart.git
cd asgart
git submodule init
git submodule update
```

You can then build ASGART by running `cargo`, the Rust build tool

```
cargo build --release
```

Once the build is finished, you will find the binaries in `target/release/asgart-*`.


# Usage

## Simple Usage

First, let us take a look at a simple example:

```
asgart seq.fasta
```

This command will look for duplications in the `seq.fasta` file, then
write them in a JSON file in the folder from which it was launched. ASGART
will probe using 20-mers, and guarantee that no duplication will
include gaps longer than 100bp in their arm-to-arm pairwise alignment.

If you wish to search reversed-complemented duplications, use the
`-R` and `-C` options, that can be combined in `-RC`. And the `-v`
option will give you more informations as the progress goes on.

```
asgart -RCv seq.fasta
```

## Input

As input(s), ASGART takes one or more FASTA files containing the
sequences within which to look for duplications. They can be either in
the FASTA (one sequence per file) or multiFASTA (multiple sequencesper
file) format.

## Output

ASGART will write its result in a JSON file in the folder
where it was launched, using the following structure:

```
{
        "strand": {
                "name":   the file(s) set by the user,
                "length": total length of the dataset,
                "map": [
                        {
                                "name":     FASTA fragment name,
                                "position": offset in the FASTA file,
                                "length":   FASTA fragment length
                        }
                ]
        },

        "settings": {
                "probe_size":             probe size used,
                "max_gap_size":           maximal gap size used,
                "min_duplication_length": minimal length for a duplicon,
                "max_cardinality":        maximal size of a family,
                "skip_masked":            were masked nucleotides skipped?,
                "trim":                   the start and end position in the dataset if it was trimmed,
        },

        "families": [            # all families
                        [        # one of these families
                            {    # a duplicon in this family
                                    "global_left_position":  position of the left arm in the input sequences,
                                    "global_right_position": position of the right arm in the input sequences,

                                    "chr_left":              fragment in the input containing the left arm,
                                    "chr_right":             fragment in the input containing the right arm,

                                    "chr_left_position":     position of the left arm relative to the start of its fragment,
                                    "chr_right_position":    position of the right arm relative to the start of its fragment,

                                    "left_length":           length of the left arm of the duplicon (bp),
                                    "right_length":          length of the right arm of the duplicon (bp),

                                    "reversed":              true if the duplication is reversed, false otherwise,
                                    "complemented":          true if the duplication is complemented, false otherwise,
                                    "identity":              the distance between the two duplicons (0.0 if not computed)
                            },
                            ...
                        ]
        ]
}
```

You can use the companion program `asgart-splice` to convert JSON files to another format.

## Options

### Functional

  - `--probe-size`/`-k` set the probing k-mers length (default: 20)

  - `--gap-size`/`-g` set the maximal gap length in a duplicon (default: 100)

  - `--min-length SIZE` specifies the minimal length (in bp) over
    which a duplication is kept in the final result and not discarded
    (default: 1000)

  - `--reverse`/`-R` look for reverse duplications

  - `--complement`/`-C` look for complemented duplications

  - `--skip-masked`/`-S` skip soft-masked zones, _i.e._ lowercased
    parts of the input files (default: no)

  - `--max-cardinality` specifies the maximal count of members in a
    duplication family (default: 500)

### Technical

  - `-h`, `--help` display an help screen

  - `-v`, `-vv`, `-vvv` increase verbosity level

  - `--out FILENAME` specifies the file in which to write the results

  - `--prefix NAME` defines a prefix to prepend to the standard output
    file name

  - `--threads COUNT` set the numbers of thread to use. Defaults to
    the number of cores abailable on the CPU

  - `--trim START END` run ASGART only on the specified area (in bp) of the
    dataset

# Plotting

ASGART comes with a plotting tool, producing a visual overview of the
duplications. Currently, four types of plots are available: chord
plots, flat plots, genome plots and Circos plots.

## Quick Start

`asgart-plot chr22.json chr22_RC.json flat`

## Arguments

`asgart-plot` takes two mandatory arguments:

1. one or more JSON-files containing results from ASGART runs;

2. the type of plot to generate.

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

  - `--no-intra` do not plot intra-fragment duplications

  - `--no-inter` do not plot inter-fragments duplications

  - `--features FILE` add an additional track containing features to
    plot alongside the duplications.

  - `--restrict-fragments A B ...` only plots fragments whose names
    are given

  - `--exclude-fragments A B ...` do not plot fragments whose names
    are given

  - `--filter-features DISTANCE` don't plot duplications that are
    farther away then `DISTANCE` bp from the features in the track.

  - `--min-thickness` set the minimal graphical width of a duplicon
    (default: 0.1)

  - `--colorize TYPE` set the method used to colorize the duplicons.
    Options are `by-type` (different colors for direct and palindromic
    duplications); `by-position` (color depends on the duplication
    position within the input file(s)); `by-fragment` (each
    duplication is colorized according to its left-most duplicons);
    `none` (all are drawn in medium grey).

### Features File Format

Features files can be provided in two formats, either in GFF3 files, or using a custom,
denser format described below.

The custom format features file format is made of a list of lines, one per feature, with
three semi-colons-separated values for each:

1. The label of the feature;
2. the start of the feature. It may either be a single integer
   representing its absolute coordinate, or be of the form
   `NAME+OFFSET`, defining a start position at `OFFSET` from the start
   of `NAME` chromosomes (from the input FASTA file);
3. The length of the feature in base pairs.

Comment lines starts with a `#`.

#### Example

```
# This is a comment line
# This is a feature named MYH14, 122358bp long, and starting at the 50,188,186th base of the chromosome 19
MYH14;19+50188186;122358
# This is a feature named Foo, starting on the 123,456,789th base of the input FASTA file and 1250bp long
Foo;123456789;1250
```

## Chord Plots

A chord plot represents duplications amongst a DNA fragment as arcs
linking point on a circle figuring a fragment. Their width is directly
proportional to the length of the duplicons they represent.

### Example

`asgart-plot human_genome.json chord --out=flat.svg --min-length 20000`

![Chord plot example](screenshots/chord.png)

## Flat Plot

Flat plots are made of two superposed horizontal bars, representing
the  concatenated fragments analyzed by ASGART, with lines linking left and
right parts of the duplicons found, their width being proportional to the
length of the duplicaton.

### Example

`asgart-plot human_Y.json flat --out=flat.svg --no-direct --no-uncomplemented --min-length 2000`

![Chord plot example](screenshots/flat.png)

## Genome Plot

Genome plots draw one bar split in four lanes per fragment. The two
leftmost lanes represente respectively the intrachromosomal direct and
palindromic duplications families, and the two rightmost respectively
the interchromosomal direct and palindromic duplications families.

### Example

`asgart-plot chr10-chrY.json genome --min-length 10000`

![Chord plot example](screenshots/genome.png)

## Circos Plots

ASGART can generate files that can be used as in input for the
[Circos](http://circos.ca/) plotting tool. Although the most important
files is arguably the `<out>.links` file (containing the duplicons to
plot), ASGART also generates minimal `<out>.conf` and
`<out>.karyotype` files, as to ensure a minimal working example to be
later expanded and/or customized according to your needs.

`asgart-plot` needs to refer to files found in the Circos distribution. Thus, the
`CIRCOS_ROOT` environment variable should be set to point at the root
of the Circos distribution. Otherwise, ASGART will generate an
`<out>.conf` file containing `{circos_root}` placeholders to be
manually replaced.

### Example

`asgart-plot human_Y.json human_Y_RC.json circos --min-length 10000`

# Change Log

_Please note that ASGART follows the [semver](https://semver.org/) versioning scheme, where an increase in the major version number reflects a non backward-compatible update._

## v2.3.0

- `asgart-cat` has been renamed to `asgart-splice`
- `asgart` does not feature multiple output formats anymore; `asgart-splice` is to
  to be used instead.

## v2.2.1

- Various minor refactoring & bug-fixes

## v2.2.0

- `asgart-concat` has been renamed to `asgart-cat`
- `asgart-cat` now offers filtering options
- `asgart-cat` now takes advantage of multi-cores CPU when possible
- `asgart-plot` now offers more filtering options
- `asgart-plot` now let the user customizes the minimal graphical
  width of a duplicon with `--min-thickness`
- `asgart-plot` now offer several algorithms to set duplicons colors
- Various bug-fixes

## v2.1.1

- Fix manifest file

## v2.1.0

- Ensure that multiple fragments in a mFASTA file are processed separately
- Add a flag to specify the minimum width of a chord
- Add filtering options
- Add tooltips to chord graphs
- Fix output files naming scheme

## v2.0.2

- Fix a bound-checking bug where the last chunk would not be processed.

## v2.0.1

- Fix a bug where a strand void of large N swaths would not be
  processed.

## v2.0

- ASGART does not differentiate anymore between strand A and strand B,
  but simply works on an arbitrarily large set of files. Thus, the
  user **SHOULD PROVIDE EACH FILE ONLY ONCE**. Moreover, it is not
  necessarily to concatenate multiple input files in a single one
  anymore. This **breaking change** should give more flexibility to
  the users and potentially simplifies pipeline design.
- The ASGART automaton has been redesigned from scratch to take into
  account interlaced SDs at nearly no cost in computation time. For
  this reason, interlaced duplication families research is now the
  only and default mode.
- ASGART will now ignore large expanses of nucleotides to ignore (Ns
  and/or masked ones) in processed strands, thus slightly improving
  performances.
- Taking advantage of these new features, the parallelization system
  has been rewritten to (i) introduce parallelism at the scale of the
  automaton; and (ii) make use of the “natural” aforementioned
  breakpoints as delimiters for chunks to process in parallel. By
  doing so, it is guaranteed (i) that no duplication families that
  would be situated between two chunks will be missed; (ii) that
  ASGART will make use of available cores even when processing less
  chunks than authorized threads.
- ASGART will now make use of the trimming feature to reduce memory
  consumption. The suffix array will be built only for the trimmed
  part, instead than for the whole input. The whole input will then be
  compared to the trimmed part, contrary to what happened in version
  1.x. Such an arrangement sacrifice some CPU power in exchange of a
  strongly reduced memory consumption when processing trimmed inputs.
  It can be used to process large sequences by trimming them in
  several consecutive subsequences, then merging the results later on.
- The JSON and GFF3 output formats have been modified to reflect the
  duplication families clustering. *Please note that they are thus
  incompatible with previous versions JSON files.*
- A new tool `asgart-concat` has been added to safely concatenate JSON
  files resulting from partial runs on the same dataset. Its intended
  use is to easily merge the results from multiple runs on the same
  dataset with different settings, e.g. direct & palindromic
  duplications or if the workload was divided in multiple sub-jobs
  using trimming.
- Plotting utilities have been modified to reflect these changes.
- The automaton will progressively grow the maximal gap size when
  extending large duplications, thus letting larger duplications arms
  be found in a less fragmented way.
- The logging system has been improved to be more detailed and more
  coherent in its way to present informations.
- Minor technical issues have been resolved: ASGART will correctly
  only use the `ID` field of FASTA files and not the subsequent
  informations; the progress bar does not glitch anymore.

## v1.5

- New, **non-retrocompatible** JSON output format containing positions
  of the duplicons both globally in the strand and relative to the
  fragment they are situated on
- `asgart-plot` can now superpose several files in a single plot
- ASGART can optionally compute the Levenshtein distance between
  duplicons
- User can set the chunking size for parallel processing (defaults to
  1,000,000bp)
- Improve output files naming
- Fix a bug in post-processing
- Fix several minor bugs in logging system
- Minor under-the-hood refactoring and improvements

## v1.4.0

- Add Jaccard distance computation to estimate identity between duplicons
- Increase font size for feature plotting

## v1.3.3

- Fix regression

## v1.3.2

- Fix arg name runtime error

## v1.3.1

- Fix erroneous GFF3 output: seq names are now corrent, no superfluous underscore and correct, relative positions instead of absolute ones.

## v1.3

- Add a new plot format, _genome_
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
