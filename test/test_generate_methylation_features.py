import shutil
import pytest
import polars as pl
import pandas as pd
from SemiBin.generate_methylation import *
import unittest
from unittest.mock import patch, MagicMock
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO



class SetupArgs:
    def __init__(self):
        self.pileup = "test/methylation_data/geobacillus-plasmids.pileup.bed"
        self.contig_fasta = "test/methylation_data/geobacillus-plasmids.assembly.fasta"
        self.bin_motifs = "test/methylation_data/bin-motifs.tsv"
        self.data = "test/methylation_data/data.csv"
        self.data_split = "test/methylation_data/data_split.csv"
        self.num_process = 1
        self.min_motif_observations_bin = 400
        self.min_motif_obs_contig = 1
        self.min_motif_methylation = 0.5
        self.min_valid_read_coverage= 1
        self.output = "test/methylation_data/test_output"

@pytest.fixture
def data():
    args = SetupArgs()
    logger = MagicMock()
    pileup, data, data_split, bin_consensus = load_data(args, logger)
    
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    assembly = read_fasta(args.contig_fasta)
    
    # Get the unique motifs
    motifs = bin_consensus\
        .select(["motif", "mod_position", "mod_type", "n_mod_bin", "n_nomod_bin"])\
        .with_columns(
            motif_mod = pl.col("motif") + "_" + pl.col("mod_position").cast(pl.String) + "_" + pl.col("mod_type"),
            n_motifs = pl.col("n_mod_bin") + pl.col("n_nomod_bin")
        )\
        .filter(pl.col("n_motifs") >= args.min_motif_observations_bin)\
        .unique(["motif_mod"])

    if len(motifs) == 0:
        logger.error(f"No motifs found")
        sys.exit(1)
    
    motif_list = [(Motif(row[0], row[1]), row[2])for row in motifs.unique(["motif", "mod_position", "mod_type"]).iter_rows()]
    contigs_in_split = data_split.select("contig").to_pandas()
    contigs_in_split = contigs_in_split["contig"].str.rsplit("_",n=1).str[0].unique()

    contigs = {}
    for c in data.get_column("contig"):
        if c in contigs_in_split:
            contigs[c] = True
        else:
            contigs[c] = False
    
    return {
        "args": args,
        "data_split": data_split,
        "data": data,
        "pileup": pileup,
        "bin_consensus": bin_consensus,
        "assembly": assembly,
        "contigs": contigs,
        "motif_list": motif_list
    }


def test_generate_methylation_features(data):
    logger = MagicMock()
    args = SetupArgs()

    data_before = pl.read_csv(args.data)
    generate_methylation_features(logger, args)

    assert os.path.exists(os.path.join(args.output, "data.csv"))
    assert os.path.exists(os.path.join(args.output, "contig_methylation.tsv"))

    data_after = pl.read_csv(os.path.join(args.output, "data.csv"))
    assert len(data_before.columns) != len(data_after.columns)
    
    os.remove(os.path.join(args.output, "data_split.csv"))
    os.remove(os.path.join(args.output, "data.csv"))
    os.remove(os.path.join(args.output, "contig_methylation.tsv"))
    os.remove(os.path.join(args.output, "contig_split.fasta"))


def test_split_contigs(data):
    args = data["args"]
    lengths = get_split_contig_lengths(data["assembly"], ["contig_2", "contig_3"])

    create_assembly_with_split_contigs(data["assembly"], lengths, os.path.join(args.output, "plasmid_split.fasta"))

    create_split_pileup(args.pileup, lengths, os.path.join(args.output, "plasmid_split.bed"))

    run_epimetheus(
        os.path.join(args.output, "plasmid_split.bed"),
        os.path.join(args.output, "plasmid_split.fasta"),
        ["GATC_a_1", "GATC_m_3"],
        min_valid_read_coverage = 1,
        threads = 1,
        output = os.path.join(args.output, "contig_split_methylation.tsv")        
    )

    expected_output = pd.read_csv(
        os.path.join("test/methylation_data/expected_contig_methylation.tsv"),
        sep = "\t"
    )

    current_output = pd.read_csv(
        os.path.join(args.output, "contig_split_methylation.tsv"),
        sep = "\t"
    )

    pd.testing.assert_frame_equal(expected_output, current_output)

    os.remove(os.path.join(args.output, "plasmid_split.bed"))
    os.remove(os.path.join(args.output, "plasmid_split.fasta"))
    os.remove(os.path.join(args.output, "contig_split_methylation.tsv"))
    shutil.rmtree(args.output)




        
class TestCreateAssemblyWithSplitContigs(unittest.TestCase):
    def setUp(self):
        # Set up a mock assembly with test contigs
        self.assembly = {
            "contig_1": SeqRecord(Seq("ATGCGTACGTAGCTAGCTAG"), id="contig_1"),
            "contig_2": SeqRecord(Seq("TGCATGCTAGCTGACTGACT"), id="contig_2")
        }
        
        # Define contigs to split
        self.split_contigs = ["contig_1", "contig_2"]

        # Output path for the test
        self.output_path = "test_split_contigs.fasta"

    def test_create_assembly_with_split_contigs(self):
        # Call the function to create split contigs
        contig_lengths = get_split_contig_lengths(self.assembly, self.split_contigs)
        create_assembly_with_split_contigs(self.assembly, contig_lengths, self.output_path)

        # Read the output file and check its contents
        with open(self.output_path, "r") as output_handle:
            split_records = list(SeqIO.parse(output_handle, "fasta"))
        
        # Expected results: each contig should be split into two parts
        expected_sequences = [
            ("contig_1_1", "ATGCGTACGT"),
            ("contig_1_2", "AGCTAGCTAG"),
            ("contig_2_1", "TGCATGCTAG"),
            ("contig_2_2", "CTGACTGACT")
        ]
        
        # Check that we have the expected number of records
        self.assertEqual(len(split_records), len(expected_sequences))

        # Check that each split record has the correct ID and sequence
        for record, (expected_id, expected_seq) in zip(split_records, expected_sequences):
            self.assertEqual(record.id, expected_id)
            self.assertEqual(str(record.seq), expected_seq)

    def tearDown(self):
        # Clean up the test output file
        if os.path.exists(self.output_path):
            os.remove(self.output_path)

class TestCheckFilesExist(unittest.TestCase):

    @patch('os.path.exists')
    def test_all_files_exist(self, mock_exists):
        # Setup the mock to return True for all paths
        mock_exists.return_value = True
        
        # args = SetupArgs('motifs_scored.txt', 'data.txt', 'data_split.txt', 'assembly.fasta', 'motif_index_dir')
        args = SetupArgs()
        # No exception should be raised if all files exist
        paths = [args.pileup, args.data, args.data_split, args.contig_fasta]
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



# def test_generate_methylation_features():
#     args = SetupArgs()
    
#     # create a mock logger
#     logger = MagicMock()
    
#     generate_methylation_features(logger, args)
    
#     assert os.path.exists(os.path.join(args.output, "data_split.csv")), "data_split.csv should be created."
#     assert os.path.exists(os.path.join(args.output, "data.csv")), "data.csv should be created."
    
#     # Cleanup
#     os.remove(os.path.join(args.output, "data_split.csv"))
#     os.remove(os.path.join(args.output, "data.csv"))
#     os.remove(os.path.join(args.output, "contig_methylation.tsv"))
#     os.remove(os.path.join(args.output, "contig_split.fasta"))
#     os.rmdir(args.output)



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
