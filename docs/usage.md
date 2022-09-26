# Usage examples

## Single-sample binning

Inputs required: `S1.fa` (contig file in FASTA format) and `S1.sorted.bam` (short reads mapped to the contigs, sorted).

[[How to generate inputs to SemiBin.](../generate)]


### Easy single binning mode

There are two options.

**1. Using a pre-trained model.** This is the fastest option and should work the best if you have metagenomes from one of our prebuilt habitats (alternatively, you can use the `global` "habitat" which combines all of them).


```bash
SemiBin single_easy_bin \
        --environment human_gut \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```

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
11. `global`

[Figure 5 in the manuscript](https://www.nature.com/articles/s41467-022-29843-y#Fig5) shows details of how well each habitat-specific model performs.

**2. Learn a new model.** Alternatively, you can learn a new model for your data.
The main disadvantage is that this approach will take a lot more time and use a lot more memory.
While using a pre-trained model should take a few minutes and use 4-6GB of RAM, training a new model may take several hours and use 40GB of RAM.

```bash
SemiBin single_easy_bin \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```

The [Supplemental Tables 5 & 6 in the SemiBin manuscript](https://www.nature.com/articles/s41467-022-29843-y#MOESM1) contain a lot more information with respect to the computational trade-offs.

If you have a lot of samples that are similar to each other while not fitting into any of our builtin trained models, you can also build your own model from a subset of them (scroll down to the "SemiBin(pretrained)" Section).

### Advanced single-sample binning workflows

The basic idea of using SemiBin with single-sample and co-assembly is:

1. generate _data.csv_ and _data_split.csv_ (used in training) for every sample,
2. train the model for every sample, and
3. bin the contigs with the model trained from the same sample.

You can run the individual steps by yourself, which can enable using compute clusters to make the binning process faster.

In particular, `single_easy_bin` includes the following steps:

`generate_cannot_links`,`generate_data_single` and `bin`; while `multi_easy_bin`
includes the following steps: `generate_cannot_links`, `generate_data_multi` and `bin`.

(1)  Generate features (`data.csv/data_split.csv` files)

```bash
SemiBin generate_sequence_features_single -i S1.fa -b S1.sorted.bam -o S1_output
```
(2) Generate cannot-link file
```bash
SemiBin generate_cannot_links -i S1.fa -o S1_output
```

Be warned that this will run `mmseqs2`, which takes a lot of time.
If you are running `mmseqs2` against the GTDB taxonomy as part of of your pipeline already, you can make SemiBin skip this step by passing it the results using the `--taxonomy-annotation-table` argument.

(3) Train a model (if desired)

```bash
SemiBin train \
    -i S1.fa \
    --data S1_output/train.csv \
    --data-split S1_output/train_split.csv \
    -c S1_output/cannot/cannot.txt \
    -o S1_output
```

This step heavily benefits from having access to a GPU.
It will attempt to auto-detect whether one is available, but you can use the `--engine` argument to specify whether SemiBin should use CPU or GPU processing.

This step can be skipped if you want to use a pretrained model.

(4) Bin
```bash
SemiBin bin \
    -i S1.fa \
    --model S1_output/model.h5 \
    --data S1_output/data.csv \
    -o S1_output
```
or, to use one of our built-in models (see above for the list of available models), you replace the `--model` argument with the `--environment` argument:

```bash
SemiBin bin \
    -i S1.fa \
    --environment human_gut \
    --data S1_output/data.csv \
    -o S1_output
```

### SemiBin(pretrain)

Another suggestion is that you can pre-train a model from part of your dataset, which can provide a balance as it is faster than training for each sample while achieving better results than a pre-trained model from another dataset (see the [manuscript](https://www.nature.com/articles/s41467-022-29843-y) for more information).

If you have `S1.fa`, `S1/data.csv`, `S1/data_split.csv`, `S1/cannot/cannot.txt`; `S2.fa`, `S2/data.csv`, `S2/data_split.csv`, `S2/cannot/cannot.txt`; `S3.fa`, `S3/data.csv`, `S3/data_split.csv`, `S3/cannot/cannot.txt`.
You can train the model from 3 samples.

```bash
SemiBin train \
    -i S1.fa S2.fa S3.fa \
    --data S1/train.csv S2/train.csv S3/train.csv \
    --data-split S1/train_split.csv S2/train_split.csv S3/train_split.csv \
    -c S1/cannot.txt s2/cannot.txt S3/cannot.txt \
    --mode several \
    -o S1_output
```

## Co-assembly binning

Input: `contig.fa` and `S1.sorted.bam`, `S2.sorted.bam`, `S3.sorted.bam`,...

### Easy co-assembly binning mode

To a large extent, co-assembly binning is just like single-sample binning.
The major difference is that when generating features, we can use multiple samples.
Unfortunately, this also means that prebuilt models cannot be used (because models depend on the number of samples, one would need to pre-train a model for each possible input sample input number).


```bash
SemiBin single_easy_bin \
    -i contig.fa \
    -b S1.sorted.bam S2.sorted.bam S3.sorted.bam \
    -o co-assembly_output
```

### Advanced co-assembly binning workflows


(1)  Generate `data.csv/data_split.csv`

```bash
SemiBin generate_sequence_features_single \
    -i contig.fa \
    -b S1.sorted.bam S2.sorted.bam S3.sorted.bam \
    -o contig_output
```

Note that we use the `generate_sequence_features_single` mode because co-assembly and single-sample modes are very similar.

(2) Generate cannot-link
```bash
SemiBin generate_cannot_links -i contig.fa -o contig_output
```
(3) Train
```bash
SemiBin train \
    -i contig.fa \
    --data contig_output/train.csv \
    --data-split contig_output/train_split.csv \
    -c contig_output/cannot/cannot.txt \
    -o contig_output
```

SemiBin will attempt to detect a GPU and fallback to CPU if none is found, but you can use the `--engine` argument to specify which one to use.
Having access to a GPU can speed up this mode.

(4) Bin

```bash
SemiBin bin \
    -i contig.fa \
    --model contig_output/model.h5 \
    --data contig_output/data.csv \
    -o output
```


## Multi-sample binning

Multi-sample binning requires more complex steps to prepare input data as well as more computation but can also result in more bins (particularly in complex habitats).
See [Figure 3b in the manuscript](https://www.nature.com/articles/s41467-022-29843-y#Fig3) for a comparison of multi sample vs. single sample.

Inputs:
- original FASTA files: `S1.fa`, `S2.fa`, `S3.fa`, `S4.fa`, and `S5.fa` (we will assume that there are 5 samples)
- combined FASTA file: `concatenated.fa` (can be generated with the `concatenate_fasta` SemiBin subcommand)
- mapped reads to the combined FASTA file: `S1.sorted.bam`, `S2.sorted.bam`, `S3.sorted.bam`, `S4.sorted.bam`, and `S5.sorted.bam`.


### Generating `concatenated.fa`

```bash
SemiBin concatenate_fasta \
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
SemiBin multi_easy_bin \
        -i concatenated.fa \
        -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam \
        -o multi_output
```

### Advanced multi-sample binning workflows workflows

As with the other modes, the `multi_easy_bin` subcommand encapsulates a series of steps that can also be run independently for more control.
They can also be parallelized in a compute cluster, for example.

(1)  Generate `data.csv/data_split.csv`

```bash
SemiBin generate_sequence_features_multi \
    -i concatenated.fa \
    -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam \
    -o output
```
(2) Generate cannot-link

You need to call `generate_cannot_links` for all the input FASTA files, independently:

```bash
for sample in S1 S2 S3 S4 S5; do
    SemiBin generate_cannot_links -i ${sample}.fa -o ${sample}_output
done
```

We used a bash for loop above, but it is equivalent to running the following:

```bash
SemiBin generate_cannot_links -i S1.fa -o S1_output
SemiBin generate_cannot_links -i S2.fa -o S2_output
SemiBin generate_cannot_links -i S3.fa -o S3_output
SemiBin generate_cannot_links -i S4.fa -o S4_output
SemiBin generate_cannot_links -i S5.fa -o S5_output
```

See the comment above about how you can bypass most of the computation if you have run `mmseqs2` to annotate your contigs against GTDB already.

(3) Train

Training is performed independently for each sample (thus, could be parallelized), but uses the _input features that account for all the sample data_:

```bash
for sample in S1 S2 S3 S4 S5 ; do
    SemiBin train \
        -i ${sample}.fa \
        --data multi_output/samples/${sample}/train.csv \
        --data-split multi_output/samples/${sample}/train_split.csv \
        -c ${sample}_output/cannot/cannot.txt \
        -o ${sample}_output
done
```

The same comments about GPU access that applied to the other modes, apply here.
Together with running `mmseqs2`, this is the more computational intensive step and may benefit from being run in parallel for each sample.

(4) Bin

```bash
for sample in S1 S2 S3 S4 S5 ; do
    SemiBin bin \
        -i ${sample}.fa \
        --model ${sample}_output/model.h5 \
        --data multi_output/samples/${sample}/data.csv \
        -o output
done
```

Each sample is binned independently.
This step is relatively fast.

