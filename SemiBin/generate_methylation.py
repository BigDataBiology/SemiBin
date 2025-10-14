#!/usr/bin/env python

import numpy as np
import multiprocessing
# from pymethylation_utils.utils import run_epimetheus
from epymetheus import epymetheus
import os
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
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
    

def process_contig_split_methylation(
    contig_name,
    contig_length,
    pileup_path,
    assembly_path,
    motifs,
    min_valid_read_coverage,
    min_valid_cov_to_diff_fraction,
):
    try:
        pileup_df = epymetheus.query_pileup_records(
            pileup_path=pileup_path,
            contigs=[contig_name],
        )
    except Exception as e:
        return None

    contig_half = contig_length // 2
    pileup_split_df = pileup_df\
        .with_columns(
            pl.when(pl.col("start") < contig_half).then(pl.col("contig").cast(pl.String) + "_1").otherwise(pl.col("contig").cast(pl.String) + "_2").alias("contig")
        )\
        .with_columns(
            pl.when(pl.col("start") < contig_half).then(pl.col("start")).otherwise(pl.col("start") - contig_half).alias("start")
        )

    contig_meth_features = epymetheus.methylation_pattern_from_dataframe(
        pileup_df=pileup_split_df,
        assembly=assembly_path,
        motifs = motifs,
        output_type=epymetheus.MethylationOutput.Median,
        threads=1,
        min_valid_read_coverage=min_valid_read_coverage,
        min_valid_cov_to_diff_fraction=min_valid_cov_to_diff_fraction
    )

    return contig_meth_features


def find_data_split_methylation_parallel(
    contigs,
    contig_lengths,
    pileup_path,
    assembly_path,
    motifs,
    threads,
    min_valid_read_coverage,
    min_valid_cov_to_diff_fraction,
):
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(threads) as pool:
        results = pool.starmap(
            process_contig_split_methylation,
            [(
                contig,
                contig_lengths[contig],
                pileup_path,
                assembly_path,
                motifs,
                min_valid_read_coverage,
                min_valid_cov_to_diff_fraction
            ) for contig in contigs]
        )

    valid_results = [df for df in results if df is not None and not df.is_empty()]
    if not valid_results:
        return pl.DataFrame()

    combined_df = pl.concat(valid_results, how = "vertical")\
        .with_columns([
            pl.col("contig").str.slice(0, pl.col("contig").str.len_chars() - 2).alias("base_contig"),
            pl.col("contig").str.slice(-1, 1).alias("split_num")
        ])\
        .sort([
            "base_contig",
            "motif",
            "mod_position",
            "mod_type",
            "split_num"
        ]).drop(["base_contig", "split_num"])
    return combined_df
    

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

        
def check_files_exist(paths=[]):
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


def sort_columns(cols):
    mod_columns = sorted([col for col in cols if "methylation_value" in col], key=lambda x: x.split("_")[-2:])
    nomod_columns = sorted([col for col in cols if "motif_present" in col], key=lambda x: x.split("_")[-2:])
    # Interleave the mod and nomod columns
    sorted_columns = [val for pair in zip(mod_columns, nomod_columns) for val in pair]
    return ["contig"] + sorted_columns  # Keep 'contig' as the first column

def create_methylation_matrix(methylation_features):
    """
    Creates a feature matrix with methylation from motifs-scored or methylation features.
    """
    # check if the methylation features have the required columns
    required_columns = ["contig", "motif", "mod_type",  "mod_position", "methylation_value", "motif_present"]
    if not all(col in methylation_features.columns for col in required_columns):
        raise ValueError(f"Missing required columns in methylation features. Required columns: {', '.join(required_columns)}")
    
    # Calculate mean methylation for each motif
    matrix = methylation_features\
        .with_columns(
            motif_mod = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String)
        )
    
    
    matrix = matrix.select(["contig", "motif_mod", "methylation_value", "motif_present"])\
        .pivot(
            index = "contig",
            on = "motif_mod",
            values = ["methylation_value", "motif_present"],
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



def check_data_file_args(logger, data, data_split, args):
    if data and data_split:
        logger.info("Using provided data and data_split files.")
    elif data or data_split:
        logger.error("Missing data or data_split path. Either both should be provided or none.")
        sys.exit(1)
    else:
        logger.info("Using default data and data_split files. Checking output directory...")
        data = os.path.join(args.output, "data.csv")
        data_split = os.path.join(args.output, "data_split.csv")
    return data, data_split
        

def generate_methylation_features(logger, contig_fasta_path, pileup_path, args, data_path = None, data_split_path = None):
    logger.info("Adding Methylation Features")
    logger.info("Loading data...")
    
    # Check for the data and data_split file
    data_path, data_split_path = check_data_file_args(logger, data_path, data_split_path, args)
        
    paths = [pileup_path, data_path, data_split_path, contig_fasta_path]
    check_files_exist(paths)
    
    # check if output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.motifs_file:
        motifs_file = pl.read_csv(args.motifs_file, separator="\t")

        # Get the unique motifs
        motifs = motifs_file\
            .select(["motif", "mod_position", "mod_type"])\
            .with_columns(
                motif_mod = pl.col("motif") + "_" + pl.col("mod_type") + "_" + pl.col("mod_position").cast(pl.String)
            )\
            .get_column("motif_mod")\
            .unique()
    elif args.motifs:
        motifs = args.motifs
    else:
        logger.error("No motifs provided. Exiting")
        sys.exit(1)

    if len(motifs) == 0:
        logger.error(f"No motifs found")
        sys.exit(1)
        
    # Load the data
    data = pl.read_csv(data_path)
    data = data\
        .rename({"": "contig"})
    data_split = pl.read_csv(data_split_path)
    data_split = data_split\
        .rename({"": "contig"})
    
    
    # Load the assembly file
    assembly = read_fasta(contig_fasta_path)

    # create splitted assembly
    contigs_to_split = data_split.select("contig").to_pandas()
    contigs_to_split = contigs_to_split["contig"].str.rsplit("_",n=1).str[0].unique()

    contig_lengths_for_splitting = get_split_contig_lengths(assembly, contigs_to_split)
    logger.info("Splitting assembly")
    create_assembly_with_split_contigs(
        assembly, contig_lengths_for_splitting , os.path.join(args.output, "contig_split.fasta")
    )

    number_of_motifs = len(motifs)
    logger.info(f"Motifs found (#{number_of_motifs}): {motifs}")

    # Run methylation utils
    logger.info("Running epimetheus for whole contigs")
    contig_methylation = epymetheus.methylation_pattern(
        pileup = pileup_path,
        assembly = contig_fasta_path,
        output = os.path.join(args.output,"contig_methylation.tsv"),
        motifs = motifs,
        threads = args.num_process,
        min_valid_read_coverage = args.min_valid_read_coverage,
        # batch_size = 1000,
        min_valid_cov_to_diff_fraction = 0.80,
        allow_assembly_pileup_mismatch = False,
        output_type=epymetheus.MethylationOutput.Median
    )

    contig_methylation.write_csv(
        os.path.join(args.output, "contig_methylation_features.tsv"),
        separator = "\t"
    )
    
    logger.info("Running epimetheus for split contigs")
    contig_split_methylation = find_data_split_methylation_parallel(
        contigs = contigs_to_split,
        contig_lengths=contig_lengths_for_splitting,
        pileup_path=pileup_path,
        assembly_path=os.path.join(args.output, "contig_split.fasta"),
        motifs =motifs,
        min_valid_read_coverage=args.min_valid_read_coverage,
        min_valid_cov_to_diff_fraction=0.80,
        threads=args.num_process
    )

    if contig_split_methylation.is_empty():
        logger.error("Failed to generate split contig methylation features. No valid results returned.")
        raise RuntimeError("No valid methylation features generated for split contigs")

    
    contig_split_methylation.write_csv(
        os.path.join(args.output, "contig_split_methylation_features.tsv"),
        separator = "\t"
    )


    contig_methylation = contig_methylation\
        .with_columns(
            pl.when(pl.col("n_motif_obs") > 0).then(1).otherwise(0).alias("motif_present")
        )

    contig_split_methylation = contig_split_methylation\
        .with_columns(
            pl.when(pl.col("n_motif_obs") > 0).then(1).otherwise(0).alias("motif_present")
        )
    
    # Methylation number is median of mean methylation. Filtering is applied for too few motif observations.
    # contig_methylation = contig_methylation.filter((pl.col("N_motif_obs") * pl.col("mean_read_cov")) >= 16)
    # contig_split_methylation = contig_split_methylation.filter((pl.col("N_motif_obs") * pl.col("mean_read_cov")) >= 16)
    contig_methylation = contig_methylation.filter((pl.col("n_motif_obs") >= args.min_motif_observations) & (pl.col("mean_read_cov") >= args.min_valid_read_coverage))
    contig_split_methylation = contig_split_methylation.filter((pl.col("n_motif_obs") >= args.min_motif_observations) & (pl.col("mean_read_cov") >= args.min_valid_read_coverage))

    print(contig_split_methylation)
    print(contig_methylation)
    
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
    print(data_methylation_matrix)
    
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

    try:
        logger.info("Writing to data and data_split files...")
        print(data.columns)
        print("Saving new data files")
        data_split.write_csv(os.path.join(args.output, "data_split.csv"), separator=",", quote_style='never') 
        data.write_csv(os.path.join(args.output, "data.csv"), separator=",", quote_style='never')
    except Exception as e:
        print(f"An error occurred while writing the output: {e}")
        sys.exit(1)
    
    
def generate_methylation_features_multi(
    logger,
    pileup_paths,
    args,
    sample_list = None,
):
    
    # Check if sample names and pileup name match.
    import copy
    pileup_dict = {}

    for p in pileup_paths:
        if not isinstance(p, str) or not os.path.isfile(p):
            raise ValueError(f"Invalid file: {p}")

        base_name = os.path.basename(p).rsplit(".", 1)[0]

        if base_name in pileup_dict:
            logger.error(f"Duplicate entry detected ({base_name}). Exiting")
            sys.exit(1)

        pileup_dict[base_name] = p

    if not sample_list:
        samples_dir = os.path.join(args.output, "samples")
        sample_list = [d for d in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, d))]

    missing_samples = set(sample_list) - set(pileup_dict.keys())
    if missing_samples:
        raise ValueError(f"Missing pileup files for samples: {', '.join(missing_samples)}")

    for sample in sample_list:
        logger.info(f"Finding methylation pattern of: {sample}")
        run_args = copy.copy(args)
        run_args.output = os.path.join(args.output, "samples", sample)
        generate_methylation_features(
            logger = logger,
            contig_fasta_path = os.path.join(args.output, "samples", f"{sample}.fa"),
            pileup_path = pileup_dict[sample],
            args = run_args,
            data_path = os.path.join(run_args.output, "data.csv"),
            data_split_path = os.path.join(run_args.output, "data_split.csv")
        )

