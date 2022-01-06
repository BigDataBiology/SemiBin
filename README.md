# SemiBin: Semi-supervised Metagenomic Binning Using Siamese Neural Networks

Command tool for metagenomic binning with semi-supervised deep learning using
information from reference genomes.

[![Test Status](https://github.com/BigDataBiology/SemiBin/actions/workflows/semibin_test.yml/badge.svg)](https://github.com/BigDataBiology/SemiBin/actions/workflows/semibin_test.yml)
[![Documentation Status](https://readthedocs.org/projects/semibin/badge/?version=latest)](https://semibin.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

_NOTE_: This tool is still in development. You are welcome to try it out and
feedback is appreciated, but expect some bugs/rapid changes until it
stabilizes. Please use [Github
issues](https://github.com/BigDataBiology/SemiBin/issues) for bug reports and
the [SemiBin users mailing-list](https://groups.google.com/g/semibin-users) for
more open-ended discussions or questions.

If you use this software in a publication please cite:

> SemiBin: Incorporating information from reference genomes with
> semi-supervised deep learning leads to better metagenomic assembled genomes
> (MAGs)
> Shaojun Pan, Chengkai Zhu, Xing-Ming Zhao, Luis Pedro Coelho
> bioRxiv 2021.08.16.456517; doi:
> [https://doi.org/10.1101/2021.08.16.456517](https://doi.org/10.1101/2021.08.16.456517)


## Install

SemiBin runs on Python 3.7-3.9.

### Install from bioconda ### 

_Note_ : If you want to use SemiBin with GPU, you need to install Pytorch with GPU support. Or `conda install -c bioconda semibin` just install Pytorch with CPU support.

```bash
conda create -n SemiBin python==3.7
conda activate SemiBin
conda install -c conda-forge -c bioconda semibin
conda install pytorch torchvision torchaudio cudatoolkit=10.2 -c pytorch-lts
```

### Install from source

You can download the source code from github and install.

Install dependence packages using conda:
[MMseqs2](https://github.com/soedinglab/MMseqs2),
[Bedtools](http://bedtools.readthedocs.org/]), [Hmmer](http://hmmer.org/),
[Fraggenescan](https://sourceforge.net/projects/fraggenescan/).

```bash
conda install -c conda-forge -c bioconda mmseqs2=13.45111
conda install -c bioconda bedtools hmmer fraggenescan==1.30
```

Once the dependencies are installed, you can install by running:

```bash
python setup.py install
```

## Examples

**NOTE**: The `SemiBin` API is a work-in-progress. The examples refer to
version `0.5`, but this may change in the near future (after the release of
version 1.0, we expect to freeze the API for [at least 5
years](https://big-data-biology.org/software/commitments/). We are very happy
to [hear any feedback on API
design](https://groups.google.com/g/semibin-users), though.

SemiBin runs on single-sample, co-assembly and multi-sample binning. 

The basic idea of using SemiBin with single-sample and co-assembly is:

(1) Generate _data.csv_ and _data_split.csv_(used in training) for every sample

(2) Training the model for every sample

(3) Binning the contigs for every contig with the model trained from the same sample

When using multi-sample binning, the basic idea is very similar, the inputs are the contigs combined from several samples and bam files from severl samples. And then we also generated _data.csv_ and _data_split.csv_, training and binning for every  sample. The only difference compared to single-sample binning is the _data.csv_ and _data_split.csv_ has the abundance information from several samples.

Considering the issue that contig annotations and model training requires significant computational time and the algorithm design of SemiBin, we proposed SemiBin(pretrain) for _single-sample binning_: 

(1) Trained a model from one sample or several samples (Or used our built-in pretrained model)

(2) Directly applied this model to other samples.

For the details and examples of every command to run SemiBin with these binning modes,  please read [read the docs](https://semibin.readthedocs.io/en/latest/usage/).

## Easy single/co-assembly binning mode

You will need the following inputs:

1. A contig file (`contig.fna` in the example below)
2. BAM files from mapping short reads to the contigs

The `single_easy_bin` command can be used in single-sample and co-assembly
binning modes (contig annotations using mmseqs with GTDB reference genome) to
get results in a single step.

```bash
SemiBin single_easy_bin -i contig.fna -b *.bam -o output --recluster
```

In this example, SemiBin will download GTDB to
`$HOME/.cache/SemiBin/mmseqs2-GTDB/GTDB`. You can change this default using the
`-r` argument.

You can set `--recluster` to use the reclustering step with single-copy genes
described in the paper, which can make results a little better (especially when the number of samples used is larger 5).

You can use `--environment` with (`human_gut`, `dog_gut`, or `ocean`) to use one of our built-in models. (**Note:** Recommended way, which will save much time for contig annotations and model training, and also get very good results) 

```bash
SemiBin single_easy_bin -i contig_S1.fna -b S1.bam -o output --environment human_gut --recluster
```

## Easy multi-samples binning mode

The `multi_easy_bin` command can be used in multi-samples binning modes (contig
annotations using mmseqs with GTDB reference genome).


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

You can get the results with one line of code. You can set `--recluster` to use
the reclustering part with single-copy genes described in the paper.

```bash
SemiBin multi_easy_bin -i contig_whole.fna -b *.bam -o output --recluster
```

## Output

The output folder will contain

1. Datasets used for training and clustering.

2. Saved semi-supervised deep learning model.

3. Output bins.

4. Some intermediate files.

For every sample, reconstructed bins are in `output_bins` directory. Using
reclustering, reconstructed bins are in `output_recluster_bins` directory.

For more details about the output, [read the
docs](https://semibin.readthedocs.io/en/latest/output/).

## Advanced workflows

You can run individual steps by yourself, which can enable using compute
clusters to make the binning process faster (especially in multi-samples
binning mode). For example, `single_easy_bin` includes the following steps:
`predict_taxonomy`,`generate_data_single` and `bin`; while `multi_easy_bin`
includes following step: `predict_taxonomy`, `generate_data_multi` and `bin`.

In advanced mode, you can also use our built-in pre-trained model in
single-sample binning mode. Here we provide pre-trained models for human gut,
dog gut and marine environment. You can just use these models for single-sample
binning and it will save much time for contig annotations and model training.

A very easy way to run SemiBin with a built-in model
(`human_gut`/`dog_gut`/`ocean` for single-sample binning):

```bash
SemiBin single_easy_bin -i contig_S1.fna -b S1.bam -o output --environment human_gut --recluster
```

Another suggestion is that you can pre-train a model from part of your dataset,
which can provide a balance as it is faster than training for each sample while
achieving better results than a pre-trained model from another dataset.

(1) Generate `data.csv/data_split.csv` for every sample

Single-sample/co-assembly binning:

```bash
SemiBin generate_data_single -i contig_S1.fna -b S1.bam -o output
```

Multi-sample binning:

```bash
SemiBin generate_data_multi -i contig_combined.fna -b S1.bam S2.bam S3.bam S4.bam S5.bam -o output -s :
```

(2) Generate cannot-link for every sample

```bash
SemiBin predict_taxonomy -i contig_S1.fna -o output
```

(3) Train a pre-trained model across several samples (For single-sample binning, Make sure the input files are corresponding)

```bash
SemiBin train -i S1.fna S2.fna S3.fna --data S1/train.csv S2/train.csv S3/train.csv --data-split S1/train_split.csv S2/train_split.csv S3/train_split.csv -c S1/cannot.txt s2/cannot.txt S3/cannot.txt -o output --mode several 
```
Or just train a model from one sample. If you are using multi-sample binning, here the contig is the contig for every sample.

Single-sample/co-assembly binning:

```bash
SemiBin train -i contig.fna --data train.csv --data-split train_split.csv -c cannot.txt -o output --mode single
```

Multi-sample binning(similar for other samples) :

```bash
SemiBin train -i contig_S1.fna --data S1/train.csv --data-split S1/train_split.csv -c S1/cannot.txt -o output --mode single --recluster
```

(4) Bin with the trained model.  If you are using multi-sample binning, here the contig is the contig for every sample and the bam files are still from several samples.

Single-sample/co-assembly binning:

```bash
SemiBin bin -i contig.fna --model model.h5 --data data.csv -o output --recluster
```

Multi-sample binning(similar for other samples):

```bash
SemiBin bin -i contig_S1.fna --model model.h5 --data S1/data.csv -o output --recluster
```

Or our built-in model(human_gut, dog_gut or ocean) (Just for single-sample binning)

```bash
SemiBin bin -i contig.fna --data data.csv -o output --environment human_gut --recluster
```

For more details on usage, including every command on how to run SemiBin with different binning modes, please [read the docs](https://semibin.readthedocs.io/en/latest/usage/).

