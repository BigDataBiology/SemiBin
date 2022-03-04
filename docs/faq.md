# Frequently Asked Questions

## Can I use another version of GTDB for annotation?

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

## Does SemiBin work with long-read data?

Technically, yes, you can apply it to long-read data and it will produce bins.
However, SemiBin is not optimized for this setting and all the benchmarking in the manuscript is performed on short-read assemblies. You may consider using SemiBin as part of a multi-algorithm approach followed by dereplication, but on its own it will be likely outperformed by methods specifically addressing long-read data (e.g., [GraphMB](https://doi.org/10.1101/2022.02.25.481923)).

How to adapt the approach of SemiBin to long-read data is part of ongoing research.

