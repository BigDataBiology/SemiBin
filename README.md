# S<sup>3</sup>N<sup>2</sup>Bin (Semi-supervised Siamese Neural Network for metagenomic binning)

[![Test Status](https://github.com/BigDataBiology/S3N2Bin/actions/workflows/s3n2bin_test.yml/badge.svg)](https://github.com/BigDataBiology/S3N2Bin/actions/workflows/s3n2bin_test.yml)
[![Documentation Status](https://readthedocs.org/projects/s3n2bin/badge/?version=latest)](https://s3n2bin.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

_NOTE_: This tool is still in development. You are welcome to try it out and
feedback is appreciated, but expect some bugs/rapid changes until it
stabilizes.

Command tool for metagenomic binning with semi-supervised deep learning using
information from reference genomes.

## Install

S<sup>3</sup>N<sup>2</sup>Bin runs on Python 3.6-3.8 and requires
[Bedtools](https://github.com/arq5x/bedtools2) to be available for calculating
depth from bam files.

### Install from source

You can download the source code from github and install with the standard

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
with GTDB reference genome). `single_easy_bin` includes the following parts: `predict_taxonomy`,`generate_data_single` and `bin`.

```bash
S3N2Bin single_easy_bin -i contig.fna -b *.bam -r /mmseqs_data/GTDB -o output
```

If you do not set the path of GTDB, S<sup>3</sup>N<sup>2</sup>Bin will download GTDB  to $HOME/.cache/S3N2Bin/mmseqs2-GTDB/GTDB.

## Easy multi-samples binning mode

The `multi_easy_bin` command can be used in
multi-samples binning modes (contig annotations using mmseqs
with GTDB reference genome). `multi_easy_bin` includes following parts: `predict_taxonomy`, `generate_data_multi` and `bin`.

You will need the following inputs.

1. A combined contig file 

2. BAM files from mapping

  For every contig, format of the name is <sample_name>:<contig_name>, : is the separator. You can use any separator you want by set `--separator` . *Note:* Make sure the sample names are unique and  the separator does not introduce confusion when splitting. For example:

```bash
<S1>:<Contig_1>
ATGCAAAA
<S1>:<Contig_2>
ATGCAAAA
<S1>:<Contig_3>
ATGCAAAA
<S2>:<Contig_1>
ATGCAAAA
<S2>:<Contig_2>
ATGCAAAA
<S3>:<Contig_1>
ATGCAAAA
<S3>:<Contig_2>
ATGCAAAA
<S3>:<Contig_3>
ATGCAAAA
```

You can get the results with one line of code. 

```bash
S3N2Bin multi_easy_bin -i contig_whole.fna -b *.bam -r /mmseqs_data/GTDB -o output -s :
```

## Advanced-bin mode

You can run every step by yourself, which will make the binning process a bit faster especially in multi-samples binning mode.

## Generate Cannot-link contrains

You can use [mmseqs2](https://github.com/soedinglab/MMseqs2) or
[CAT](https://github.com/dutilh/CAT) (or other contig annotation tools) to get
taxonomic classifications of contigs. 

If you want to use mmseqs(default in S<sup>3</sup>N<sup>2</sup>Bin), you can use subcommand `predict_taxonomy` to generate the cannot-link file. If you want to use CAT, you can run CAT first and  use the script
`script/concatenate.py` to generate the cannot-link file(contig1, contig2) that
can be used in S<sup>3</sup>N<sup>2</sup>Bin.

#### mmseqs2

```bash
S3N2Bin predict_taxonomy -i contig.fna -r /mmseqs_data/GTDB -o output
```

#### CAT

```bash
python script/concatenate.py -i CAT.out -c contig.fna -s sample-name -o output --CAT
```

## Generating data for training and clustering(data.csv;data_split.csv)

### Single/co-assembly binning

```bash
S3N2Bin generate_data_single -i contig.fna -b *.bam -o output
```

### Multi-samples binning

```bash
S3N2Bin generate_data_multi -i contig_whole.fna -b *.bam -o output -s :
```

## Binning(training and clustering)

If you run S<sup>3</sup>N<sup>2</sup>Bin on a GPU server, S<sup>3</sup>N<sup>2</sup>Bin will run on GPU automatically.

```bash
S3N2Bin bin -i contig.fna --data data.csv --data-split data_split.csv -c cannot.txt -o output
```

For more details of the usage, please  [read the docs](https://s3n2bin.readthedocs.io/en/latest/usage/). 

## Output

The output folder will contain

1. Datasets used for training and clustering.

2. Saved semi-supervised deep learning model.

3. Output bins.

4. Some intermediate files.

For every sample, reconstructed bins are in `output_recluster_bins` directory.

For more details about the output, [read the docs](https://s3n2bin.readthedocs.io/en/latest/output/). 

