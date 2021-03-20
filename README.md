# S³N²Bin (Semi-supervised Siamese Neural Network for metagenomic binning)

[![Test Status](https://github.com/BigDataBiology/S3N2Bin/actions/workflows/s3n2bin_test.yml/badge.svg)](https://github.com/BigDataBiology/S3N2Bin/actions/workflows/s3n2bin_test.yml)
[![Documentation Status](https://readthedocs.org/projects/s3n2bin/badge/?version=latest)](https://s3n2bin.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

_NOTE_: This tool is still in development. You are welcome to try it out and
feedback is appreciated, but expect some bugs/rapid changes until it
stabilizes. Please use [Github
issues](https://github.com/BigDataBiology/S3N2Bin/issues) for bug reports and
the [Discussions](https://github.com/BigDataBiology/S3N2Bin/discussions) for
more open-ended discussions/questions.

Command tool for metagenomic binning with semi-supervised deep learning using
information from reference genomes.

## Install

S<sup>3</sup>N<sup>2</sup>Bin runs on Python 3.6-3.8.

### Install from source

You can download the source code from github and install.

Install dependence packages using conda: [Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/),  [Fraggenescan](https://sourceforge.net/projects/fraggenescan/) and [cmake](https://cmake.org/).

```bash
conda install -c bioconda bedtools hmmer fraggenescan
```
```bash
conda install -c anaconda cmake=3.19.6
```

```bash
python setup.py install
```

## Examples

## Easy single/co-assembly binning mode

You will need the following inputs:

1. A contig file (`contig.fna` in the example below)
2. BAM files from mapping

You can get the results with one line of code. The `single_easy_bin` command can be used in
single-sample and co-assembly binning modes (contig annotations using mmseqs
with GTDB reference genome). `single_easy_bin` includes the following steps:
`predict_taxonomy`,`generate_data_single` and `bin`.

```bash
S3N2Bin single_easy_bin -i contig.fna -b *.bam -o output
```

In this example, S³N²Bin will download GTDB to
`$HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB`. You can change this default using the
`-r` argument.

## Easy multi-samples binning mode

The `multi_easy_bin` command can be used in
multi-samples binning modes (contig annotations using mmseqs
with GTDB reference genome). `multi_easy_bin` includes following step:
`predict_taxonomy`, `generate_data_multi` and `bin`.

You will need the following inputs.

1. A combined contig file

2. BAM files from mapping

For every contig, format of the name is `<sample_name>:<contig_name>`, where
`:` is the default separator (it can be changed with the `--separator`
argument). *Note:* Make sure the sample names are unique and  the separator
does not introduce confusion when splitting. For example:

```bash
>S1:Contig_1
AGATAATAAAGATAATAATA
>S1:Contig_2
CGAATTTATCTCAAGAACAAGAAAA
>S1:Contig_3
AAAAAGAGAAAATTCAGAATTAGCCAATAAAATA
>S2:Contig_1
AATGATATAATACTTAATA
>S2:Contig_2
AAAATATTAAAGAAATAATGAAAGAAA
>S3:Contig_1
ATAAAGACGATAAAATAATAAAAGCCAAATCCGACAAAGAAAGAACGG
>S3:Contig_2
AATATTTTAGAGAAAGACATAAACAATAAGAAAAGTATT
>S3:Contig_3
CAAATACGAATGATTCTTTATTAGATTATCTTAATAAGAATATC
```

You can get the results with one line of code.

```bash
S3N2Bin multi_easy_bin -i contig_whole.fna -b *.bam -o output
```

## Advanced-bin mode

You can run individual steps by yourself, which can enable using compute
clusters to make the binning process faster (especially in multi-samples
binning mode).

For more details on usage, including information on how to run individual steps
separately, [read the docs](https://s3n2bin.readthedocs.io/en/latest/usage/).

## Output

The output folder will contain

1. Datasets used for training and clustering.

2. Saved semi-supervised deep learning model.

3. Output bins.

4. Some intermediate files.

For every sample, reconstructed bins are in `output_recluster_bins` directory.

For more details about the output, [read the docs](https://s3n2bin.readthedocs.io/en/latest/output/).

