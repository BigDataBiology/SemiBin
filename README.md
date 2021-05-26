# SemiBin (Semi-supervised Siamese Neural Network for metagenomic binning)

[![Test Status](https://github.com/BigDataBiology/SemiBin/actions/workflows/semibin_test.yml/badge.svg)](https://github.com/BigDataBiology/SemiBin/actions/workflows/semibin_test.yml)
[![Documentation Status](https://readthedocs.org/projects/semibin/badge/?version=latest)](https://semibin.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

_NOTE_: This tool is still in development. You are welcome to try it out and
feedback is appreciated, but expect some bugs/rapid changes until it
stabilizes. Please use [Github
issues](https://github.com/BigDataBiology/SemiBin/issues) for bug reports and
the [Discussions](https://github.com/BigDataBiology/SemiBin/discussions) for
more open-ended discussions/questions.

Command tool for metagenomic binning with semi-supervised deep learning using
information from reference genomes.

## Install

SemiBin runs on Python 3.6-3.9.

### Install from source

You can download the source code from github and install.

Install dependence packages using conda: [MMseqs2](https://github.com/soedinglab/MMseqs2), [Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/),  [Fraggenescan](https://sourceforge.net/projects/fraggenescan/) and [cmake](https://cmake.org/).

```bash
conda install -c conda-forge -c bioconda mmseqs2=13.45111
conda install -c bioconda bedtools hmmer fraggenescan
conda install -c anaconda cmake=3.19.6
```

Once the dependencies are installed, you can install by running:

```bash
python setup.py install
```

## Examples

## Advanced-bin mode

You can run individual steps by yourself, which can enable using compute
clusters to make the binning process faster (especially in multi-samples
binning mode). 

In advanced-bin mode, you can also use our built-in pre-trained model in single-sample binning mode. Here we provide pre-trained models for human gut, dog gut and tara environment. You can just use these models for single-sample binning and it will save much time for contig annotations and model training. 

In our experiments, we found that training for every sample then binning would get the best results, but it need much time. Using our provided trained model is a good option and it can also get very good results and significantly perform better than Metabat2. 

A very easy way to run SemiBin with a built-in model([human_gut/dog_gut/ocean]):

```bash
SemiBin single_easy_bin -i contig.fna -b *.bam -o output --environment human_gut
```

Another suggestion for running with pre-trained model is that you can get a pre-trained model from your dataset, which is a better way that we think can save time and get good results at the same time. For example, you can subsample several samples(i.e. 10 or 15) as training samples. Then you can train a good pre-trained model from these sample and then you can use this pre-trained model for binning. 

(1) Generate data.csv/data_split.csv for every sample

```bash
SemiBin generate_data_single -i contig.fna -b *.bam -o output
```

(2) Generate cannot-link for every sample

```bash
SemiBin predict_taxonomy -i contig.fna -o output
```

(3) Train a pre-trained model across several samples

```bash
SemiBin train -i *.fna --data *.csv --data-split *.csv -c cannot*.txt -o output --mode several
```

(4) Bin with the trained model

```bash
SemiBin bin -i contig.fna --model model.h5 --data data.csv -o output 
```

For more details on usage, including information on how to run individual steps
separately, [read the docs](https://semibin.readthedocs.io/en/latest/usage/).

## Easy single/co-assembly binning mode

You will need the following inputs:

1. A contig file (`contig.fna` in the example below)
2. BAM files from mapping

You can get the results with one line of code. The `single_easy_bin` command can be used in
single-sample and co-assembly binning modes (contig annotations using mmseqs
with GTDB reference genome). `single_easy_bin` includes the following steps:
`predict_taxonomy`,`generate_data_single` and `bin`.

```bash
SemiBin single_easy_bin -i contig.fna -b *.bam -o output
```

In this example, SemiBin will download GTDB to
`$HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB`. You can change this default using the
`-r` argument. You can set `--recluster` to use the reclustering part with single-copy genes described in the paper, which can make results a little better.

You can use `--environment` with(human_gut, dog_gut and ocean) to use our built-in model.

```bash
SemiBin single_easy_bin -i contig.fna -b *.bam -o output --environment human_gut
```

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

You can get the results with one line of code. You can set `--recluster` to use the reclustering part with single-copy genes described in the paper.

```bash
SemiBin multi_easy_bin -i contig_whole.fna -b *.bam -o output
```

## Output

The output folder will contain

1. Datasets used for training and clustering.

2. Saved semi-supervised deep learning model.

3. Output bins.

4. Some intermediate files.

For every sample, reconstructed bins are in `output_bins` directory. Using reclustering, reconstructed bins are in `output_recluster_bins` directory.

For more details about the output, [read the docs](https://semibin.readthedocs.io/en/latest/output/).

