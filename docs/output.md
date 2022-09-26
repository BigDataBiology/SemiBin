# Outputs of SemiBin

## Single sample/co-assembly binning

* `output_recluster_bins`: directory of all reconstructed bins after reclustering.
* `output_bins`: directory of all reconstructed bins before reclustering.
* `model.h5`: saved semi-supervised deep learning model. 
* `data.csv/data_split.csv`: data used in the training of deep learning model.
* `*_data_cov.csv/*_data_split_cov.csv`: coverage data generated from depth file.
* `cannot/cannot.txt`: cannot-link file used in the training.
* `recluster_bins_info.tsv`: table with basic information on each bin (name, total number of basepairs, number of contigs, N50, and L50; more columns may be added in the future)

## Multi-samples binning

* `bins`: Reconstructed bins from all samples.
* `samples/*.fasta`: Contig fasta file for every sample from the input whole_contig.fna.
* `samples/*_data_cov.csv`: same in single sample/coassembly binning.
* `samples/{sample-name}/`: directory of the output of SemiBin for every sample(same as that in Single sample/coassembly binning). 

