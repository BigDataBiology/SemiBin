# S<sup>3</sup>N<sup>2</sup>Bin (Semi-supervised Siamese Neural Network for metagenomic binning)

[![License: MIT](https://anaconda.org/bioconda/gmgc-mapper/badges/license.svg)](https://anaconda.org/bioconda/gmgc-mapper)

Command tool for metagenomic binning with semi-supervised deep leaning using additional information from reference genomes.

## Install

$S^3N^2Bin$ runs on Python 3.6-3.8 and requires [Bedtools] (https://github.com/arq5x/bedtools2) to be available for calculating depth from bam files.

### Install from source

You can download the source code from github and install with the standard

```bash
python setup.py install
```

## Generate Cannot-link contrains

  You can use [CAT](https://github.com/dutilh/CAT)(or other contig annotation tools) to get the taxonomic classifications of contigs. Then you can use the script `script/concatenate.py` to generate the cannot-link file(contig_1<TAB>contig_2 ) that can be used in $S^3N^2Bin$. 

```bash
python script/concatenate.py -i CAT.out -c contig.fna -s sample-name -o output
```

## Examples

### Single sample/coassembly binning

```bash
S3N2Bin -i contig.fna -b *.bam -c cannot-link.txt -o output 
```

### Multiple samples binning(Must set -s parameter)

```bash
S3N2Bin -i whole_contig.fna -b *.bam -c *.txt -s C -o output
```

#### Multiple samples binning pipeline

(1) Concatenate all contigs from all samples together. Make sure that names of samples are unique and id for every contig is <sample_name><separator><contig_id>. Concatenated contig format is:

```bash
<S1>O<Contig_1>
ATGCAAAA
<S1>O<Contig_2>
ATGCAAAA
<S1>O<Contig_3>
ATGCAAAA
<S2>O<Contig_1>
ATGCAAAA
<S2>O<Contig_2>
ATGCAAAA
<S3>O<Contig_1>
ATGCAAAA
<S3>O<Contig_2>
ATGCAAAA
<S3>O<Contig_3>
ATGCAAAA
```

(2) Map reads to the concatenated contig file to get the bam files.

(3) Generate cannot-link files for every sample. The name of the cannot-link file is sample_name.txt. Make sure the sample name here is same to that in step(1).

(4) Run $S^3N^2Bin$

## Output

The output folder will contain

1. Datasets used for training and clustering.

2. Saved semi-supervised deep learning model.

3. Output bins.

4. Some intermediate files. 

For more details, read the docs.

