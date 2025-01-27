#!/usr/bin/env python

import numpy as np
from multiprocessing import get_context
import multiprocessing
from pymethylation_utils.utils import run_epimetheus
import os
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


# os.environ['POLARS_MAX_THREADS'] = str(args.num_process)
import polars as pl


def read_fasta(path):
    # Check if the file exists
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")    
    
    # Check if the file has a valid FASTA extension
    valid_extensions = ['.fasta', '.fa', '.fna', '.gz']
    if not any(path.endswith(ext) for ext in valid_extensions):
        raise ValueError(f"Unsupported file extension. Please provide a FASTA file with one of the following extensions: {', '.join(valid_extensions)}")
    
    # Check if the file is a gzipped FASTA file
    if path.endswith('.gz'):
        with gzip.open(path, "rt") as handle:  # "rt" mode for reading as text
            contigs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
                    
    else:
        # Read a regular (uncompressed) FASTA file
        with open(path, "r") as handle:
            contigs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return contigs

def get_split_contig_lengths(assembly, split_contigs):
    contig_lengths = {}
    for c in split_contigs:
        contig = assembly[c]
        contig_len = len(contig.seq)
        contig_lengths[contig.id] = contig_len

    return contig_lengths
    

def create_assembly_with_split_contigs(assembly, contig_lengths, output):
    split_records = []
    for c in contig_lengths.keys():
        contig = assembly[c]

        contig_half = contig_lengths[c] // 2

        s1 = contig.seq[:contig_half]
        s2 = contig.seq[contig_half:]

        r1 = SeqRecord(
            s1,
            contig.id + "_1",
            description = ""
        )
        r2 = SeqRecord(
            s2,
            contig.id + "_2",
            description = ""
        )

        split_records.append(r1)
        split_records.append(r2)
    with open(output, "w") as output_handle:
        SeqIO.write(split_records, output_handle, "fasta")

def create_split_pileup(
    pileup_path: str,
    contig_lengths: dict[str, int],
    output_path: str
):
    """
    Reads a large TSV/CSV file line by line, splitting contigs in half
    and adjusting 'start' positions accordingly, then writes out the result.

    :param pileup_path: Path to the pileup file (tab-delimited).
    :param contig_lengths: Dictionary {contig_name: length}.
    :param output_path: Path to write the processed file.
    """

    with open(pileup_path, "r") as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            # Split into the 18 expected columns (assuming no header).
            (
                contig,
                start,
                end,
                mod_type,
                score,
                strand,
                start2,
                end2,
                color,
                N_valid_cov,
                percent_modified,
                N_modified,
                N_canonical,
                N_other_mod,
                N_delete,
                N_fail,
                N_diff,
                N_nocall,
            ) = line.strip().split("\t")

            # Filter for relevant contigs.
            if contig not in contig_lengths:
                continue

            c_len = contig_lengths[contig]

            # Convert start to an integer before comparisons.
            start_val = int(start)
            half_length = c_len // 2

            # Decide if row belongs to the first half or the second half.
            if start_val < half_length:
                contig = f"{contig}_1"
            else:
                contig = f"{contig}_2"
                start_val -= half_length

            # Write out the updated line.
            out_cols = [
                contig,
                str(start_val),
                end,
                mod_type,
                score,
                strand,
                start2,
                end2,
                color,
                N_valid_cov,
                percent_modified,
                N_modified,
                N_canonical,
                N_other_mod,
                N_delete,
                N_fail,
                N_diff,
                N_nocall,
            ]
            f_out.write("\t".join(out_cols) + "\n")
        
def check_files_exist(paths=[], directories=[]):
    """
    Checks if the given files and directories exist.
    
    Parameters:
    paths (list): List of file paths to check.
    directories (list): List of directory paths to check.
    
    Raises:
    FileNotFoundError: If any of the specified files or directories do not exist.
    """
    for f in paths:
        if not os.path.exists(f):
            raise FileNotFoundError(f"The file {f} does not exist.")
    
    for d in directories:
        if not os.path.exists(d):
            raise FileNotFoundError(f"The directory {d} does not exist.")


def load_data(args, logger):
    """Loads data"""
    pileup = pl.scan_csv(args.pileup, separator = "\t", has_header = False, new_columns = [
                             "contig", "start", "end", "mod_type", "score", "strand", "start2", "end2", "color", "N_valid_cov", "percent_modified", "N_modified", "N_canonical", "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall"
                         ])
    bin_consensus = pl.read_csv(args.bin_motifs, separator="\t")
    
    data = pl.read_csv(args.data)
    data = data\
        .rename({"": "contig"})
    data_split = pl.read_csv(args.data_split)
    data_split = data_split\
        .rename({"": "contig"})

    logger.info("Data loaded successfully.")
    return pileup, data, data_split, bin_consensus


def sort_columns(cols):
    mod_columns = sorted([col for col in cols if "median" in col], key=lambda x: x.split("_")[-2:])
    nomod_columns = sorted([col for col in cols if "motif_present" in col], key=lambda x: x.split("_")[-2:])
    # Interleave the mod and nomod columns
    sorted_columns = [val for pair in zip(mod_columns, nomod_columns) for val in pair]
    return ["contig"] + sorted_columns  # Keep 'contig' as the first column

def create_methylation_matrix(methylation_features):
    """
    Creates a feature matrix with methylation from motifs-scored or methylation features.
    """
    # check if the methylation features have the required columns
    required_columns = ["contig", "motif", "mod_type",  "mod_position", "median", "motif_present"]
    if not all(col in methylation_features.columns for col in required_columns):
        raise ValueError(f"Missing required columns in methylation features. Required columns: {', '.join(required_columns)}")
    
    # Calculate mean methylation for each motif
    matrix = methylation_features\
        .with_columns(
            motif_mod = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String)
        )
    
    
    matrix = matrix.select(["contig", "motif_mod", "median", "motif_present"])\
        .pivot(
            index = "contig",
            on = "motif_mod",
            values = ["median", "motif_present"],
            aggregate_function = None,
            maintain_order = True
        )\
        .rename(
            lambda column_name: column_name.replace("motif_mod_", "")
        )\
        .fill_null(0)


    new_columns=sort_columns(matrix.columns)
    matrix = matrix.select(new_columns)
    
    return matrix



def check_data_file_args(logger, args):
    if args.data and args.data_split:
        logger.info("Using provided data and data_split files.")
    elif args.data or args.data_split:
        logger.error("Missing data or data_split path. Either both should be provided or none.")
        sys.exit(1)
    else:
        logger.info("Using default data and data_split files. Checking output directory.")
        args.data = os.path.join(args.output, "data.csv")
        args.data_split = os.path.join(args.output, "data_split.csv")
    return args
        

def generate_methylation_features(logger, args):
    logger.info("Adding Methylation Features")    
    logger.info("Loading data...")
    
    # Check for the data and data_split file
    args = check_data_file_args(logger, args)
        
        
    paths = [args.pileup, args.data, args.data_split, args.contig_fasta, args.bin_motifs]
    directories = []

    check_files_exist(paths, directories)
    
    # check if output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Load the data
    logger.info("Loading methylation data...")
    lf_pileup, data, data_split, bin_consensus = load_data(args, logger)
    
    
    # Load the assembly file
    assembly = read_fasta(args.contig_fasta)

    # create splitted assembly
    contigs_to_split = data_split.select("contig").to_pandas()
    contigs_to_split = contigs_to_split["contig"].str.rsplit("_",n=1).str[0].unique()

    contig_lengths_for_splitting = get_split_contig_lengths(assembly, contigs_to_split)
    create_assembly_with_split_contigs(
        assembly, contig_lengths_for_splitting , os.path.join(args.output, "contig_split.fasta")
    )
    create_split_pileup(args.pileup, contig_lengths_for_splitting, os.path.join(args.output, "pileup_split.bed"))

    
        
    # Get the unique motifs
    motifs = bin_consensus\
        .select(["motif", "mod_position", "mod_type", "n_mod_bin", "n_nomod_bin"])\
        .with_columns(
            motif_mod = pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.String) ,
            n_motifs = pl.col("n_mod_bin") + pl.col("n_nomod_bin")
        )\
        .filter(pl.col("n_motifs") >= args.min_motif_observations_bin)\
        .get_column("motif_mod")\
        .unique()

    if len(motifs) == 0:
        logger.error(f"No motifs found")
        sys.exit(1)
    
        
    number_of_motifs = len(motifs)
    logger.info(f"Motifs found (#{number_of_motifs}): {motifs}")

    # Run methylation utils
    code = run_epimetheus(
        args.pileup,
        args.contig_fasta,
        motifs,
        args.num_process,
        args.min_valid_read_coverage,
        os.path.join(args.output,"contig_methylation.tsv")
    )
    if code != 0:
        logger.error("Error running methylation_utils for all contigs")
        sys.exit(1)

    
    code = run_epimetheus(
        os.path.join(args.output, "pileup_split.bed"),
        os.path.join(args.output, "contig_split.fasta"),
        motifs,
        args.num_process,
        args.min_valid_read_coverage,
        os.path.join(args.output,"contig_split_methylation.tsv")
    )
    if code != 0:
        logger.error("Error running methylation_utils for split contigs")
        sys.exit(1)

    schema = {
        'contig': pl.String(),
        'motif': pl.String(),
        'mod_type': pl.String(),
        'mod_position': pl.Int8(),
        'median': pl.Float64(),
        'mean_read_cov': pl.Float64(),
        'N_motif_obs': pl.Int32(),
        'motif_occurences_total': pl.Int32(),
    }
    contig_methylation = pl.read_csv(os.path.join(args.output, "contig_methylation.tsv"), separator = "\t", schema = schema)\
        .with_columns(
            pl.when(pl.col("N_motif_obs") > 0).then(1).otherwise(0).alias("motif_present")
        )
    contig_split_methylation = pl.read_csv(os.path.join(args.output, "contig_split_methylation.tsv"), separator = "\t", schema = schema)\
        .with_columns(
            pl.when(pl.col("N_motif_obs") > 0).then(1).otherwise(0).alias("motif_present")
        )
    

    # Methylation number is median of mean methylation. Filtering is applied for too few motif observations.
    contig_methylation = contig_methylation.filter(pl.col("N_motif_obs") >= args.min_motif_obs_contig)
    contig_split_methylation = contig_split_methylation.filter(pl.col("N_motif_obs") >= args.min_motif_obs_contig)
    
    data_split_methylation_matrix = create_methylation_matrix(
        methylation_features = contig_split_methylation
    )
    
    data_split = data_split\
        .join(
            data_split_methylation_matrix,
            on = "contig",
            how = "left"
        )\
        .rename({"contig": ''})\
        .fill_nan(0.0)\
        .fill_null(0.0)
        
    data_methylation_matrix = create_methylation_matrix(
        methylation_features=contig_methylation
    ).select(data_split_methylation_matrix.columns)
    
    data = data\
        .join(
            data_methylation_matrix,
            on = "contig",
            how = "left"
        )\
        .rename({"contig": ''})\
        .fill_nan(0.0)\
        .fill_null(0.0)
 
    assert data_split_methylation_matrix.columns == data_methylation_matrix.columns, "methylation columns does not match between data_split_methylation_matrix and data_methylation_matrix"

    data_split_cols = data_split.columns
    data_cols = data.columns

    data_cols_filtered = [col for col in data_cols if not "mean" in col]
    data_cols_filtered = [col for col in data_cols_filtered if not "var" in col]
    assert data_split_cols == data_cols_filtered, "data.csv and data_split.csv columns does not match after methylation addition"
    
    try:
        logger.info("Writing to data and data_split files...")
        data_split.write_csv(os.path.join(args.output, "data_split.csv"), separator=",", quote_style='never') 
        data.write_csv(os.path.join(args.output, "data.csv"), separator=",", quote_style='never')
    except Exception as e:
        print(f"An error occurred while writing the output: {e}")
        sys.exit(1)
    
    
    
