import pytest
from SemiBin.generate_methylation import *
import unittest
from unittest.mock import patch, MagicMock
import os


class SetupArgs:
    def __init__(self):
        self.bed = "test/methylation_data/pileup.bed"
        self.contig_fasta = "test/methylation_data/assembly.fasta"
        self.motifs_scored = "test/methylation_data/motifs-scored.tsv"
        self.bin_motifs = "test/methylation_data/bin-motifs.tsv"
        self.data = "test/methylation_data/data.csv"
        self.data_split = "test/methylation_data/data_split.csv"
        self.motif_index_dir = "test/methylation_data/motif-positions"
        self.min_fraction = 0.6
        self.num_process = 1
        self.motif_occurence_cutoff = 0.9
        self.min_motif_observations = 1
        self.output = "test_output"

@pytest.fixture
def data():
    args = SetupArgs()
    logger = MagicMock()
    motifs_scored, data, data_split, bin_consensus = load_data(args, logger)
    
    contigs = get_contigs(data_split)
    
    assembly = read_fasta(args.contig_fasta, contigs)
    
    # Get lenths of the contigs
    contig_lengths = {contig: len(sequence) for contig, sequence in assembly.items()}
    
    return {
        "args": args,
        "data_split": data_split,
        "data": data,
        "bin_consensus": bin_consensus,
        "contig_fasta": assembly,
        "motifs_scored": motifs_scored,
        "contigs": contigs,
        "contig_lengths": contig_lengths
    }


def test_get_motifs(data):
    """
    Test get_motifs function at different cutoffs.
    """
    bin_consensus = data["bin_consensus"]
    data = data["motifs_scored"]
    
    motifs_all = get_motifs(data, bin_consensus, occurence_cutoff=0)
    motifs_all_len = len([motif for motifs in motifs_all.values() for motif in motifs])
    
    motifs_90 = get_motifs(data, bin_consensus, occurence_cutoff=0.9)
    motifs_90_len = len([motif for motifs in motifs_90.values() for motif in motifs])
    
    assert motifs_all_len > motifs_90_len, "All motifs should be greater than 90% cutoff."
    assert motifs_all_len == 17, "All motifs should be 17."
    



def test_data_split_methylation_parallel(data):
    contig_lengths = data["contig_lengths"]
    motifs_scored = data["motifs_scored"]
    bin_consensus = data["bin_consensus"]
    args = SetupArgs()
    
    motifs = get_motifs(motifs_scored, bin_consensus, 0.9)
    
    contig_split_methylation = data_split_methylation_parallel(contig_lengths, motifs, args.motif_index_dir)
    
    contig_split_methylation = contig_split_methylation\
        .with_columns(
            motif = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String)
        )
        
        
    motifs_scored = motifs_scored\
        .with_columns(
            motif = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String)
        )\
        .filter(pl.col("contig") == "contig_10")
    
    
    motif = "RGATCY_a-1"
    
    split_n_mod = contig_split_methylation\
        .filter(pl.col("motif") == motif)\
        .get_column("n_mod").to_list()
        
    scored_n_mod = motifs_scored\
        .filter(pl.col("motif") == motif)\
        .get_column("n_mod").to_list()
    
    
    assert sum(split_n_mod) == sum(scored_n_mod)
    
    split_n_nomod = contig_split_methylation\
        .filter(pl.col("motif") == motif)\
        .get_column("n_nomod").to_list()
    
    scored_n_nomod = motifs_scored\
        .filter(pl.col("motif") == motif)\
        .get_column("n_nomod").to_list()
    
    assert sum(split_n_nomod) == sum(scored_n_nomod)
    



class TestCheckFilesExist(unittest.TestCase):

    @patch('os.path.exists')
    def test_all_files_exist(self, mock_exists):
        # Setup the mock to return True for all paths
        mock_exists.return_value = True
        
        # args = SetupArgs('motifs_scored.txt', 'data.txt', 'data_split.txt', 'assembly.fasta', 'motif_index_dir')
        args = SetupArgs()
        # No exception should be raised if all files exist
        paths = [args.motifs_scored, args.data, args.data_split, args.contig_fasta]
        try:
            check_files_exist(paths)
        except FileNotFoundError:
            self.fail("FileNotFoundError raised unexpectedly!")

    @patch('os.path.exists')
    def test_file_does_not_exist(self, mock_exists):
        # Setup the mock to return False when checking for the first missing file
        def side_effect(arg):
            return arg != 'data.txt'
        
        mock_exists.side_effect = side_effect
        
        # args = SetupArgs('motifs_scored.txt', 'data.txt', 'data_split.txt', 'assembly.fasta', 'motif_index_dir')
        args = SetupArgs()
        paths = ["data.txt"]
        with self.assertRaises(FileNotFoundError) as context:
            check_files_exist(paths)
        
        # Check if the error message is correct
        self.assertIn('The file data.txt does not exist.', str(context.exception))

    @patch('os.path.exists')
    def test_directory_does_not_exist(self, mock_exists):
        # Assume all files exist but the directory does not
        def side_effect(arg):
            if arg == 'motif_index_dir':
                return False
            return True

        mock_exists.side_effect = side_effect
        
        dirs = ['motif_index_dir']
        with self.assertRaises(FileNotFoundError) as context:
            check_files_exist(directories = dirs)
        
        # Check if the correct exception for the directory is raised
        self.assertIn('The directory motif_index_dir does not exist.', str(context.exception))



def test_generate_methylation_features():
    args = SetupArgs()
    
    # create a mock logger
    logger = MagicMock()
    
    generate_methylation_features(logger, args)
    
    assert os.path.exists(os.path.join(args.output, "data_split.csv")), "data_split.csv should be created."
    assert os.path.exists(os.path.join(args.output, "data.csv")), "data.csv should be created."
    
    # Cleanup
    os.remove(os.path.join(args.output, "data_split.csv"))
    os.remove(os.path.join(args.output, "data.csv"))
    os.rmdir(args.output)



class TestCheckFilesAndLog(unittest.TestCase):

    @patch('sys.exit')
    def test_check_data_files_missing_data_split(self, mock_exit):
        args = SetupArgs()
        args.data_split = None
        logger = MagicMock()
        
        check_data_file_args(logger, args)
        
        # Ensure sys.exit(1) was called
        mock_exit.assert_called_once_with(1)
        
        # Ensure the correct error message was logged
        logger.error.assert_called_with("Missing data or data_split path. Either both should be provided or none.")

    @patch('sys.exit')
    def test_check_data_files_missing_data(self, mock_exit):
        args = SetupArgs()
        args.data = ""
        logger = MagicMock()
        
        check_data_file_args(logger, args)
        
        # Ensure sys.exit(1) was called
        mock_exit.assert_called_once_with(1)
        
        # Ensure the correct error message was logged
        logger.error.assert_called_with("Missing data or data_split path. Either both should be provided or none.")


    def test_check_data_files_default_files_missing(self):
        args = SetupArgs()
        args.data = None
        args.data_split = None
        logger = MagicMock()
        
        check_data_file_args(logger, args)
        
        # Ensure the correct error message was logged
        logger.info.assert_called_with("Using default data and data_split files. Checking output directory.")




if __name__ == '__main__':
    unittest.main()