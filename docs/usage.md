# Usage

## Single-sample binning

Inputs required: `S1.fa` (contig file in FASTA format) and `S1.sorted.bam` (short reads mapped to the contigs, sorted).
[How to generate inputs to SemiBin.](generate.html)


### Easy single binning mode

There are two options.

**1. Using a pre-trained model.** This is the fastest option and should work the best if you have metagenomes from one of our prebuilt habitats (alternatively, you can use the `global` "habitat" which combines all of them).
Supported habitats are (names should be self-explanatory, except `global` which is a generic model):

1.  `human_gut`
2.  `dog_gut`
3.  `ocean`
4.  `soil`
5.  `cat_gut`
6.  `human_oral`
6.  `mouse_gut`
7.  `pig_gut`
8.  `built_environment`
9.  `wastewater`
10. `global`

```bash
SemiBin single_easy_bin \
        --environment human_gut \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```

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

### Advanced workflows

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

This will run `mmseqs`, which takes a lot of time.

(3) Train
```bash
SemiBin train -i S1.fa --data S1_output/train.csv --data-split S1_output/train_split.csv -c S1_output/cannot/cannot.txt -o S1_output --mode single
```

(4) Bin
```bash
SemiBin bin -i S1.fa --model S1_output/model.h5 --data S1_output/data.csv -o output
```
or with our built-in models (see above for the list of available models)
```bash
SemiBin bin -i S1.fa --data S1_output/data.csv -o output --environment human_gut
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
    -o output
```

## Co-assembly binning

Input: `contig.fa` and `S1.sorted.bam`, `S2.sorted.bam`, `S3.sorted.bam`,...

### Easy single binning mode

This is similar to what is done for single-sampling binning, but pre-built models cannot be used.

```bash
SemiBin single_easy_bin -i contig.fa -b S1.sorted.bam S2.sorted.bam S3.sorted.bam -o output
```
### Advanced workflows

(1)  Generate `data.csv/data_split.csv` 
```bash
SemiBin generate_sequence_features_single -i contig.fa -b S1.sorted.bam S2.sorted.bam S3.sorted.bam -o contig_output
```
(2) Generate cannot-link 
```bash
SemiBin generate_cannot_links -i contig.fa -o contig_output
```
(3) Train
```bash
SemiBin train -i contig.fa --data contig_output/train.csv --data-split contig_output/train_split.csv -c contig_output/cannot/cannot.txt -o contig_output --mode single
```
(4) Bin

```bash
SemiBin bin -i contig.fa --model contig_output/model.h5 --data contig_output/data.csv -o output
```


## Multi-sample binning

Input: 
original fasta: S1.fa S2.fa S3.fa S4.fa S5.fa 
combined: combined.fa and S1.sorted.bam, S2.sorted.bam, S3.sorted.bam, S4.bam, S5.sorted.bam

The format of combined.fa: for every contig, format of the name is `<sample_name>:<contig_name>`, where
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
CAAAT
```
### Easy multi binning mode
```bash
SemiBin multi_easy_bin -i combined.fa -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam -o multi_output
```

### Advanced workflows

(1)  Generate `data.csv/data_split.csv` 
```bash
SemiBin generate_sequence_features_multi -i combined.fa -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam -o output -s :
```
(2) Generate cannot-link 
```bash
SemiBin generate_cannot_links -i S1.fa -o S1_output
```
```bash
SemiBin generate_cannot_links -i S2.fa -o S2_output
```
```bash
SemiBin generate_cannot_links -i S3.fa -o S3_output
```
```bash
SemiBin generate_cannot_links -i S4.fa -o S4_output
```
```bash
SemiBin generate_cannot_links -i S5.fa -o S5_output
```
(3) Train
```bash
SemiBin train -i S1.fa --data multi_output/samples/S1/train.csv --data-split multi_output/samples/S1/train_split.csv -c S1_output/cannot/cannot.txt -o S1_output --mode single
```
```bash
SemiBin train -i S2.fa --data multi_output/samples/S2/train.csv --data-split multi_output/samples/S2/train_split.csv -c S2_output/cannot/cannot.txt -o S2_output --mode single
```
```bash
SemiBin train -i S3.fa --data multi_output/samples/S3/train.csv --data-split multi_output/samples/S3/train_split.csv -c S3_output/cannot/cannot.txt -o S3_output --mode single
```
```bash
SemiBin train -i S4.fa --data multi_output/samples/S4/train.csv --data-split multi_output/samples/S4/train_split.csv -c S4_output/cannot/cannot.txt -o S4_output --mode single
```
```bash
SemiBin train -i S5.fa --data multi_output/samples/S5/train.csv --data-split multi_output/samples/S5/train_split.csv -c S5_output/cannot/cannot.txt -o S5_output --mode single
```
(4) Bin
```bash
SemiBin bin -i S1.fa --model S1_output/model.h5 --data multi_output/samples/S1/data.csv -o output 
```
```bash
SemiBin bin -i S2.fa --model S2_output/model.h5 --data multi_output/samples/S2/data.csv -o output
```
```bash
SemiBin bin -i S3.fa --model S3_output/model.h5 --data multi_output/samples/S3/data.csv -o output 
```
```bash
SemiBin bin -i S4.fa --model S4_output/model.h5 --data multi_output/samples/S4/data.csv -o output 
```
```bash
SemiBin bin -i S5.fa --model S5_output/model.h5 --data multi_output/samples/S5/data.csv -o output 
```

