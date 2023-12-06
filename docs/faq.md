# Frequently Asked Questions

## Does SemiBin work with long-read data?

Starting in _version 1.4_, **yes**!

While the original version of SemiBin could be applied to long-read data, it was suboptimal as the data looks different enough from short-read data (even post assembly).
Starting with version 1.4, SemiBin supports an alternative clustering procedure which gives much better results when applied to long-reads assemblies (including hybrid assemblies).
Use the flag `--sequencing-data=long_reads` when binning.

See the [SemiBin2 preprint](https://doi.org/10.1101/2023.01.09.523201) for a description and benchmarking of the long-read algorithm.

## What should I do if I have hybrid data (short- and long-reads)?

From SemiBin's point-of-view, you should generally treat this using the
long-reads pipeline (`--sequencing-type=long_read`).

## Does SemiBin work for eukaryotic genomes?

Technically, yes, you can apply it to eukaryotic data and it will produce bins.
However, SemiBin is not optimized for this setting and all the benchmarking in the manuscript is performed on prokaryotic data.
The long-read sequencing algorithm is particularly discouraged as it relies on prokaryotic single-copy genes.
Similarly, for short-reads, reclustering should be turned off (using `--no-recluster`) as it relies on the same set of genes.
You may consider using SemiBin as part of a multi-algorithm approach followed by dereplication, but we do not have the data to recommend its use on its own.

We are very keen to test SemiBin for these data and ask that, if you have eukaryotic metagenomics data, you feel free to get in touch ([shaojun@big-data-biology.org](mailto:shaojun@big-data-biology.org) or [luispedro@big-data-biology.org](mailto:luispedro@big-data-biology.org)).

## Can I use another version of GTDB for annotation?

**Note**: this is only relevant for the now-deprecated _SemiBin1_ pipeline

Yes. There are two approaches:

1. Download an mmseqs-formatted GTDB (the command `mmseqs databases GTDB GTDB
   tmp` will download the latest version). Then, point SemiBin to this database
   using the `--reference-db-data-dir` option.
2. Precompute the contig annotations with mmseqs using any version of GTDB and
   pass the contig annotation table to SemiBin using the
   `taxonomy-annotation-table` option. Do note that the tool expects an mmseqs
   formatted file and is likely to produce nonsensical results if a different
   format is provided.

The second approach is more complex but can make sense as part of a larger
pipeline where taxonomic annotation of contigs is performed for multiple
reasons (not only for the benefit of SemiBin).

