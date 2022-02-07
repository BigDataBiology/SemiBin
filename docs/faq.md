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

