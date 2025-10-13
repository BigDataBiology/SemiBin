# Running semi-supervised mode

Semi-supervised mode was the original training mode in [SemiBin1](https://www.nature.com/articles/s41467-022-29843-y). It is not recommended anymore, as it is slower and uses more memory than the self-supervised mode introudced in [SemiBin2](https://academic.oup.com/bioinformatics/article/39/Supplement_1/i21/7210480).

However, it is still available in SemiBin2 for backwards compatibility reasons. This page provides a guide on how to run SemiBin2 in semi-supervised mode, assuming you are already familiar with the basics of SemiBin2 (see the [usage](usage) page for a general overview).

## Easy modes

The easiest way to run SemiBin2 is to use the `single_easy_bin` subcommand, which will take care of all the steps for you.
If you want to train a model using the semi-supervised mode, you can use the `--semi-supervised` flag:


```bash
SemiBin2 single_easy_bin \
        --semi-supervised \
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
```


## Advanced single-sample binning workflows

Generally, this follows the same steps as in the [self-supervised mode](self-supervised), but with a few differences.

1. `generate_data_single` or `generate_data_multi` (as in the self-supervised mode)
2. `generate_cannot_links` (this is the extra step, compared to the self-supervised mode, _see below_)
3. `train_semi` (replacing `train_self` from the self-supervised mode)
4. `bin_short` or `bin_long` (as in the self-supervised mode)

### Generate cannot-link file (single or co-assembly modes)

```bash
SemiBin2 generate_cannot_links -i S1.fa -o S1_output
```

Be warned that this will run `mmseqs2`, which takes a lot of time.
If you are running `mmseqs2` against the GTDB taxonomy as part of of your pipeline already, you can make SemiBin2 skip this step by passing it the results using the `--taxonomy-annotation-table` argument.

### Generate cannot-link file (multi-sample mode)


You need to call `generate_cannot_links` for all the input FASTA files, independently:

```bash
for sample in S1 S2 S3 S4 S5; do
    SemiBin2 generate_cannot_links -i ${sample}.fa -o ${sample}_output
done
```

We used a bash for loop above, but it is equivalent to running the following:

```bash
SemiBin2 generate_cannot_links -i S1.fa -o S1_output
SemiBin2 generate_cannot_links -i S2.fa -o S2_output
SemiBin2 generate_cannot_links -i S3.fa -o S3_output
SemiBin2 generate_cannot_links -i S4.fa -o S4_output
SemiBin2 generate_cannot_links -i S5.fa -o S5_output
```

See the comment above about how you can bypass most of the computation if you have run `mmseqs2` to annotate your contigs against GTDB already.
