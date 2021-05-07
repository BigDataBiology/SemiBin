# Output of SemiBin

## Single sample/co-assembly binning

* `output_recluster_bins`: directory of all reconstructed bins after reclustering.
* `output_bins`: directory of all reconstructed bins before reclustering.
* `model.h5`: saved semi-supervised deep learning model. 
* `data.csv/data_split.csv`: data used in the training of deep learning model.
* `*_depth.txt`: depth file generated from bam file using `bedtools genomecov`.
* `*_data_cov.csv/*_data_split_cov.csv`: coverage data generated from depth file.
* `.faa/.ffn/.gff/.out/.hmmout/.hmmout.out/`: intermediate files when estimating number of single-copy marker genes.
* `.seed`: seed contig by estimating single-copy marker genes.
* `cannot/cannot.txt`: cannot-link file used in the training.

## Multi-samples binning

* `bins`: Reconstructed bins from all samples.
* `samples/*.fasta`: Contig fasta file for every sample from the input whole_contig.fna.
* `samples/*_depth.txt`: same as that in Single sample/coassembly binning.
* `samples/*_data_cov.csv`: same as that in Single sample/coassembly binning.
* `samples/{sample-name}/`: directory of the output of SemiBin for every sample(same as that in Single sample/coassembly binning). 

