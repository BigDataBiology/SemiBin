# Usage examples

## Binning modes

SemiBin2 supports three different binning modes, with different tradeoffs.

### Single-sample binning

Single sample binning means that each sample is assembled and binned independently.

This mode allows for parallel binning of samples and avoids cross-sample chimeras, but it does not use co-abundance information across samples.

Using a prebuilt model means that SemiBin2 can return results in a few minutes.

### Co-assembly binning

Co-assembly binning means samples are co-assembled first (as if the pool of samples were a single sample) and then bins are constructed from this pool of co-assembled contigs.

This mode can potentially generate better contigs (especially from species that are at a low abundance in any individual sample) and uses co-abundance information which can lead to better bins.
On the other hand, co-assembly can lead to inter-sample chimeric contigs and binning based on co-assembly does not retain sample-specific variation.

It is most appropriate when the samples are very similar and can be expected to contain overlapping sets of organisms (_e.g._, a time-series from the same habitat).

### Multi-sample binning

With multi-sample binning, multiple samples are assembled and binned individually, but _information from multiple samples is used together_.
This mode can use co-abundance information and retain sample-specific variation at the same time.
As we document in the [SemiBin1](https://www.nature.com/articles/s41467-022-29843-y) and [SemiBin2](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i21/7210480) manuscripts, this mode often returns the highest-number of bins (particularly for complex environments, such as soil).

However, it has increased computational costs.
In particular, prebuilt models cannot be used.

This mode is implemented by concatenating the contigs assembled from the individual samples together and then mapping reads from each sample to this concatenated database.
Concatenating the inputs can be done with the `concatenate_fasta` subcommand.

## Single-sample binning

Inputs required: `S1.fa` (contig file in FASTA format) and `S1.sorted.bam` (short reads mapped to the contigs, sorted).

[[How to generate inputs to SemiBin](../generate)]

### Easy single binning mode

There are two options.

**1. Using a pre-trained model.** This is the fastest option and should work the best if you have metagenomes from one of our prebuilt habitats (alternatively, you can use the `global` "habitat" which combines all of them).

```bash
SemiBin2 single_easy_bin \
        --environment human_gut \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```

Binning assemblies from long reads:

```bash
SemiBin2 single_easy_bin \
        --environment human_gut \
        --sequencing-type long_read \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```

📝 For hybrid assemblies (both short and long reads), treat them as long-read assemblies and use the `--sequencing-type long_read` flag.

Supported habitats are (names should be self-explanatory, except `global` which is a generic model):

1.  `human_gut`
2.  `dog_gut`
3.  `ocean`
4.  `soil`
5.  `cat_gut`
6.  `human_oral`
7.  `mouse_gut`
8.  `pig_gut`
9.  `built_environment`
10. `wastewater`
11. `chicken_caecum`
12. `global`

[Figure 5 in the SemiBin1 manuscript](https://www.nature.com/articles/s41467-022-29843-y#Fig5) shows details of how well each habitat-specific model performs (except for the `chicken_caecum` model which was contributed after publication by [Florian Plaza Oñate](https://scholar.google.com/citations?user=-gE5y_4AAAAJ) and is available since version 1.2).

**2a. Learn a new model (self-supervised mode).** You can also learn a new model for your data.
It will take a bit of time, but may produce better results:

```bash
SemiBin2 single_easy_bin \
        --self-supervised \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```

The [Supplemental Tables 5 & 6 in the SemiBin1 manuscript](https://www.nature.com/articles/s41467-022-29843-y#MOESM1) contain a lot more information with respect to the computational trade-offs.

If you have a lot of samples that are similar to each other while not fitting into any of our builtin trained models, you can also build your own model from a subset of them (see [[training a SemiBin model](training)])

### Advanced single-sample binning workflows

The basic pipeline using SemiBin2 for either single-sample and co-assembly modes:

1. generate _data.csv_ and _data_split.csv_ (used in training) for every sample,
2. train the model for every sample, and
3. bin the contigs with the model trained from the same sample.

You can run the individual steps by yourself, which can enable using compute clusters to make the binning process faster.

In particular, `single_easy_bin` includes the following steps:

1. `generate_sequence_features_single`
2. `train_self`
3. `bin_short` or `bin_long`

`multi_easy_bin` includes
1. `generate_sequence_features_multi`
2. `train_self` (if needed)
3. `bin_short` or `bin_long`

(1)  Generate features (`data.csv/data_split.csv` files)

```bash
SemiBin2 generate_sequence_features_single -i S1.fa -b S1.sorted.bam -o S1_output
```
(3) Train a model (if desired)

```bash
SemiBin2 train_self \
    --data S1_output/data.csv \
    --data-split S1_output/data_split.csv \
    -o S1_output
```

This step heavily benefits from having access to a GPU.
It will attempt to auto-detect whether one is available, but you can use the `--engine` argument to specify whether SemiBin2 should use CPU or GPU processing.

This step can be skipped if you want to use a pretrained model.

(4) Bin
```bash
SemiBin2 bin_short \
    -i S1.fa \
    --model S1_output/model.pt \
    --data S1_output/data.csv \
    -o S1_output
```
or
```bash
SemiBin2 bin_long \
    -i S1.fa \
    --model S1_output/model.pt \
    --data S1_output/data.csv \
    -o S1_output
```

or, to use one of our built-in models (see above for the list of available models), you replace the `--model` argument with the `--environment` argument:

```bash
SemiBin2 bin_short \
    -i S1.fa \
    --environment human_gut \
    --data S1_output/data.csv \
    -o S1_output
```
or
```bash
SemiBin2 bin_long \
    -i S1.fa \
    --environment human_gut \
    --data S1_output/data.csv \
    -o S1_output
```

### SemiBin2(pretrain)

Another suggestion is that you can pre-train a model from part of your dataset, which can provide a balance as it is faster than training for each sample while achieving better results than a pre-trained model from another dataset (see the [SemiBin1 manuscript](https://www.nature.com/articles/s41467-022-29843-y) for more information).

If you have `S1.fa`, `S1/data.csv`, `S1/data_split.csv`; `S2.fa`, `S2/data.csv`, `S2/data_split.csv`; `S3.fa`, `S3/data.csv`, `S3/data_split.csv`.
You can train the model from 3 samples.

```bash
SemiBin2 train_self \
    --train-from-many \
    -i S1.fa S2.fa S3.fa \
    --data S1/data.csv S2/data.csv S3/data.csv \
    --data-split S1/data_split.csv S2/data_split.csv S3/data_split.csv \
    -o S1_output
```

## Co-assembly binning

Input: `contig.fa` and `S1.sorted.bam`, `S2.sorted.bam`, `S3.sorted.bam`,...

### Easy co-assembly binning mode

To a large extent, co-assembly binning is just like single-sample binning.
The major difference is that when generating features, we can use multiple samples.
Unfortunately, this also means that prebuilt models cannot be used (because models depend on the number of samples, one would need to pre-train a model for each possible input sample input number).


```bash
SemiBin2 single_easy_bin \
    -i contig.fa \
    -b S1.sorted.bam S2.sorted.bam S3.sorted.bam \
    -o co-assembly_output
```

### Advanced co-assembly binning workflows


(1)  Generate `data.csv/data_split.csv`

```bash
SemiBin2 generate_sequence_features_single \
    -i contig.fa \
    -b S1.sorted.bam S2.sorted.bam S3.sorted.bam \
    -o contig_output
```

Note that we use the `generate_sequence_features_single` mode because co-assembly and single-sample modes are very similar.

(2) Train
```bash
SemiBin2 train_self \
    --data contig_output/data.csv \
    --data-split contig_output/data_split.csv \
    -o contig_output
```

SemiBin2 will attempt to detect a GPU and fallback to CPU if none is found, but you can use the `--engine` argument to specify which one to use.
Having access to a GPU can speed up this mode.

(3) Bin

```bash
SemiBin2 bin_short \
    -i contig.fa \
    --model contig_output/model.pt \
    --data contig_output/data.csv \
    -o output
```


## Multi-sample binning

Multi-sample binning requires more complex steps to prepare input data as well as more computation but can also result in more bins (particularly in complex habitats).
See [Figure 3b in the SemiBin1 manuscript](https://www.nature.com/articles/s41467-022-29843-y#Fig3) for a comparison of multi sample vs. single sample.

Inputs:
- original FASTA files: `S1.fa`, `S2.fa`, `S3.fa`, `S4.fa`, and `S5.fa` (we will assume that there are 5 samples)
- combined FASTA file: `concatenated.fa.gz` (can be generated with the `concatenate_fasta` SemiBin2 subcommand)
- mapped reads to the combined FASTA file: `S1.sorted.bam`, `S2.sorted.bam`, `S3.sorted.bam`, `S4.sorted.bam`, and `S5.sorted.bam`.


### Generating `concatenated.fa`

```bash
SemiBin2 concatenate_fasta \
    --input-fasta S1.fa S2.fa S3.fa S4.fa S5.fa \
    --output output
```

This will produce the file `output/concatenated.fa`

**Technical note on the format of `concatenated.fa`**: every contig is renamed to the name `<sample_name>:<original_contig_name>`, where `:` is the default separator (it can be changed with the `--separator` argument, which _must then be passed to all the commands that use it_).
Using the `concatenate_fasta` subcammand will make sure that sample names are unique and the separator does not introduce confusion when splitting (that is, that the separator is not already used in the contig or sample names).
Otherwise, you can also prepare the file yourself.
For example:

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
CAAAT
```

Afterwards, you should map each sample separately to the concatenated FASTA file to produce the respective `sorted.bam` file.

### Easy multi binning mode

```bash
SemiBin2 multi_easy_bin \
        -i concatenated.fa \
        -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam \
        -o multi_output
```

or for long reads:

```bash
SemiBin2 multi_easy_bin \
        --sequencing-type long_read \
        -i concatenated.fa \
        -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam \
        -o multi_output
```

📝 For hybrid assemblies (both short and long reads), treat them as long-read assemblies and use the `--sequencing-type long_read` flag.

### Advanced multi-sample binning workflows workflows

As with the other modes, the `multi_easy_bin` subcommand encapsulates a series of steps that can also be run independently for more control.
They can also be parallelized in a compute cluster, for example.

(1)  Generate `data.csv/data_split.csv`

```bash
SemiBin2 generate_sequence_features_multi \
    -i concatenated.fa.gz \
    -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam \
    -o output
```
(2) Train

Training is performed independently for each sample (thus, could be parallelized), but uses the _input features that account for all the sample data_:

```bash
for sample in S1 S2 S3 S4 S5 ; do
    SemiBin2 train_self \
        --data multi_output/samples/${sample}/data.csv \
        --data-split multi_output/samples/${sample}/data_split.csv \
        --output ${sample}_output
done
```

The same comments about GPU access that applied to the other modes, apply here.

(3) Bin

There are two subcommands, depending on whether you want to use the binning mode for short reads (`bin_short`) of for long reads (`bin_long`), so that this would be either

```bash
for sample in S1 S2 S3 S4 S5 ; do
    SemiBin2 bin_short \
        -i ${sample}.fa \
        --model ${sample}_output/model.pt \
        --data multi_output/samples/${sample}/data.csv \
        -o output
done
```
or
```bash
for sample in S1 S2 S3 S4 S5 ; do
    SemiBin2 bin_long \
        -i ${sample}.fa \
        --model ${sample}_output/model.pt \
        --data multi_output/samples/${sample}/data.csv \
        -o output
done
```

Each sample is binned independently.
This step is relatively fast.


## Running SemiBin in semi-supervised mode

**Note⚠️**: This is generally not needed as semi-supervised mode is not recommended anymore!

See the [semi-supervised mode](semi-supervised) page for more information.

## Running SemiBin with strobealign-aemb

This has its own [dedicated page](aemb).
