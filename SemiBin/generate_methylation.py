#!/usr/bin/env python

import numpy as np
from multiprocessing import get_context
import multiprocessing
from nanomotif.parallel import update_progress_bar
import os
import sys
import gzip
from Bio import SeqIO
import cProfile
import pstats
import argparse


# os.environ['POLARS_MAX_THREADS'] = '1'
import polars as pl

# IUPAC codes dictionary for sequence pattern matching
iupac_dict = {
    "A": "A", "T": "T", "C": "C", "G": "G",
    "R": "[AG]", "Y": "[CT]", 
    "S": "[GC]", "W": "[AT]", 
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ATCG]"
}


def read_fasta(path, contigs):
    # Check if the file exists
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")    
    
    # Check if the file has a valid FASTA extension
    valid_extensions = ['.fasta', '.fa', '.fna', '.gz']
    if not any(path.endswith(ext) for ext in valid_extensions):
        raise ValueError(f"Unsupported file extension. Please provide a FASTA file with one of the following extensions: {', '.join(valid_extensions)}")

    # Create a set from the contigs argument for fast lookup
    contigs_set = set(contigs)
    found_contigs = {}
    
    
    # Check if the file is a gzipped FASTA file
    if path.endswith('.gz'):
        with gzip.open(path, "rt") as handle:  # "rt" mode for reading as text
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in contigs_set:
                    found_contigs[record.id] = str(record.seq)
                    # Check if all needed contigs have been found
                    if set(found_contigs.keys()) == contigs_set:
                        break
    else:
        # Read a regular (uncompressed) FASTA file
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in contigs_set:
                    found_contigs[record.id] = str(record.seq)
                    # Check if all needed contigs have been found
                    if set(found_contigs.keys()) == contigs_set:
                        break
                    
    return found_contigs


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


def load_data(args):
    """Loads data"""
    try:
        motifs_scored = pl.read_csv(args.motifs_scored, separator="\t")
        # Remove beta initial 1 from n_nomod
        motifs_scored = motifs_scored\
            .with_columns(
                n_nomod = pl.col("n_nomod") - 1
            )
        
        bin_consensus = pl.read_csv(args.bin_motifs, separator="\t")
        
        data = pl.read_csv(args.data)
        data = data\
            .rename({"": "contig"})
        data_split = pl.read_csv(args.data_split)
        data_split = data_split\
            .rename({"": "contig"})

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)  # Exit the program with an error code

    print("All files loaded successfully.")
    return motifs_scored, data, data_split, bin_consensus


def get_contigs(data_split):
    """
    Takes the data split and returns a list of unique contigs.
    """
    contigs_in_split = data_split["contig"].str.split("_").map_elements(lambda x: "_".join(x[:-1]), return_dtype = pl.String).to_list()
    contigs_in_split = list(set(contigs_in_split))  # remove duplicates
    contigs_in_split.sort()  # Sort the list in place
    return contigs_in_split

def identify_ambiguous_motifs(motifs_scored, min_motif_observations = 8, ambiguous_interval = [0.05, 0.15], ambiguous_motif_percentage_cutoff = 0.4):
    # Remove motifs where the majority of the mean methylation of motifs is in the range of 0.05-0.4
    #
    contig_motif_mean_density = motifs_scored \
        .filter(pl.col("n_motifs") >= min_motif_observations) \
        .select(["contig", "motif_mod", "mod_position", "n_mod", "n_motifs"]) \
        .with_columns(
            (pl.col("n_mod") / pl.col("n_motifs")).alias("mean"),
        )\
        .filter(pl.col("mean") > 0.05)\
        .with_columns(
            ((pl.col("mean") > ambiguous_interval[0]) & (pl.col("mean") < ambiguous_interval[1])).alias("is_ambiguous")
        )
    
    # Identify ambigous motifs in sample
    motif_ambiguity = contig_motif_mean_density.group_by(["motif_mod"]) \
        .agg(
            pl.col("is_ambiguous").sum().alias("total_ambiguous"),
            pl.col("is_ambiguous").count().alias("n_contigs_with_motif")
        ) \
        .with_columns((pl.col("total_ambiguous") / pl.col("n_contigs_with_motif")).alias("percentage_ambiguous"))
    
    ambiguous_motifs = motif_ambiguity.filter(pl.col("percentage_ambiguous") >= ambiguous_motif_percentage_cutoff)\
        .get_column("motif_mod")
    
    return ambiguous_motifs


def get_motifs(motifs_scored, bin_consensus, occurence_cutoff=0.9, min_motif_observations = 8, ambiguous_interval = [0.05, 0.15], ambiguous_motif_percentage_cutoff = 0.4):
    """Extracts and returns unique motifs for each contig."""
    motifs_in_bin_consensus = bin_consensus\
        .select(["motif", "mod_position", "mod_type", "n_mod_bin", "n_nomod_bin"])\
        .with_columns(
            motif_mod = pl.col("motif") + "_" + pl.col("mod_position").cast(pl.String) + "_" + pl.col("mod_type"),
            n_motifs = pl.col("n_mod_bin") + pl.col("n_nomod_bin")
        )\
        .filter(pl.col("n_motifs") >= 700)\
        .get_column("motif_mod")\
        .unique()
    
    # filter motifs based on occurence.
    motifs_scored = motifs_scored\
        .with_columns(
            n_motifs = pl.col("n_mod") + pl.col("n_nomod"),
            motif_mod = pl.col("motif") + "_" + pl.col("mod_position").cast(pl.String) + "_" + pl.col("mod_type")
        )\
        .filter(pl.col("motif_mod").is_in(motifs_in_bin_consensus))
    

    # Total contigs in motifs_scored
    total_contigs_in_motifs = motifs_scored.unique(subset=["contig"]).shape[0]
    
    # TODO is n_motifs filter > 1 relevant after removing unambiguous motifs?
    # .filter(~pl.col("motif_mod").is_in(ambiguous_motifs))\
    
    motif_occurences_in_contigs = motifs_scored\
        .filter(pl.col("n_motifs") >= min_motif_observations)\
        .group_by(["motif", "mod_position", "mod_type"])\
        .agg(
            pl.count("contig").alias("motif_distinct_contigs")
        )\
        .with_columns(
            motif_occurences_per = pl.col("motif_distinct_contigs") / total_contigs_in_motifs
        ).sort("motif_distinct_contigs", descending = True)
        
    motifs = motif_occurences_in_contigs\
        .filter(pl.col("motif_occurences_per") >= occurence_cutoff)\
        .select(["motif", "mod_position", "mod_type"]).unique()
    
    motifs_dict = {}
    
    for mod_type in motifs.get_column("mod_type").unique().to_list():
        motif_types = motifs\
            .filter(pl.col("mod_type") == mod_type)\
            .select(["motif", "mod_position", "mod_type"])\
            .unique()\
            .to_dict()
            
            
        motifs_dict[mod_type] = [f"{motif}_{position}_{mod_type}" for motif, position, mod_type in zip(motif_types["motif"], motif_types["mod_position"], motif_types["mod_type"])]
    
    return motifs_dict


def find_data_split_methylation_pattern(contig_lengths, motifs, motif_index_dir):
    """
    For each contig in the data split, find the methylation pattern on each contig half.
    """
    
    methylation = pl.DataFrame()
    
    for contig in contig_lengths.keys():
        length = contig_lengths[contig]
        for mod_type in motifs.keys():
            motif_list = motifs[mod_type]
            index_file = os.path.join(motif_index_dir, f"{contig}_{mod_type}_motifs_positions.npz")
            if not os.path.exists(index_file):
                print(f"Index file not found: {index_file}")
                continue
            else:
                try:
                    with np.load(index_file, allow_pickle=True) as data:
                        data = {key: data[key] for key in data.files if key in motif_list}
                        
                        for motif in data.keys():
                            motif_data = data[motif].item()  # If it's a pickled object
                            index_meth_forward = motif_data["index_meth_fwd"]
                            index_nonmeth_forward = motif_data["index_nonmeth_fwd"]
                            index_meth_reverse = motif_data["index_meth_rev"]
                            index_nonmeth_reverse = motif_data["index_nonmeth_rev"]
        
                            n_mod = [
                                len(index_meth_forward[index_meth_forward < (length / 2)]) + len(index_meth_reverse[index_meth_reverse < (length / 2)]),
                                len(index_meth_forward[index_meth_forward >= (length / 2)]) + len(index_meth_reverse[index_meth_reverse >= (length / 2)])
                            ]
                            
                            n_nomod = [
                                len(index_nonmeth_forward[index_nonmeth_forward < (length / 2)]) + len(index_nonmeth_reverse[index_nonmeth_reverse < (length / 2)]),
                                len(index_nonmeth_forward[index_nonmeth_forward >= (length / 2)]) + len(index_nonmeth_reverse[index_nonmeth_reverse >= (length / 2)])
                            ]
                            
                            motif_str = motif.split("_")[0]
                            
                            methylation_tmp = pl.DataFrame({
                                    "contig": [f"{contig}_1", f"{contig}_2"],
                                    "motif": [motif_str, motif_str],
                                    "mod_type": [mod_type, mod_type],
                                    "mod_position": [motif.split("_")[-2], motif.split("_")[-2]],
                                    "n_mod": n_mod,
                                    "n_nomod": n_nomod
                                })
                            
                            methylation = pl.concat([methylation, methylation_tmp])
                except Exception as e:
                    print(f"Error reading index file: {e}")
                    continue
    
    return methylation


def calculate_contig_methylation_pattern(contig, contig_length, motifs, mod_type, index_file):
    
    methylation = pl.DataFrame()
    motif_list = motifs[mod_type]
    
    try:
        with np.load(index_file, allow_pickle=True) as data:
            data = {key: data[key] for key in data.files if key in motif_list}
            
            for motif in data.keys():
                motif_data = data[motif].item()
                index_meth_forward = motif_data["index_meth_fwd"]
                index_nonmeth_forward = motif_data["index_nonmeth_fwd"]
                index_meth_reverse = motif_data["index_meth_rev"]
                index_nonmeth_reverse = motif_data["index_nonmeth_rev"]

                n_mod = [
                    len(index_meth_forward[index_meth_forward < (contig_length / 2)]) + len(index_meth_reverse[index_meth_reverse < (contig_length / 2)]),
                    len(index_meth_forward[index_meth_forward >= (contig_length / 2)]) + len(index_meth_reverse[index_meth_reverse >= (contig_length / 2)])
                ]
                
                n_nomod = [
                    len(index_nonmeth_forward[index_nonmeth_forward < (contig_length / 2)]) + len(index_nonmeth_reverse[index_nonmeth_reverse < (contig_length / 2)]),
                    len(index_nonmeth_forward[index_nonmeth_forward >= (contig_length / 2)]) + len(index_nonmeth_reverse[index_nonmeth_reverse >= (contig_length / 2)])
                ]
                
                motif_str = motif.split("_")[0]
                
                methylation_tmp = pl.DataFrame({
                        "contig": [f"{contig}_1", f"{contig}_2"],
                        "motif": [motif_str, motif_str],
                        "mod_type": [mod_type, mod_type],
                        "mod_position": [motif.split("_")[-2], motif.split("_")[-2]],
                        "n_mod": n_mod,
                        "n_nomod": n_nomod
                    })
                
                methylation = pl.concat([methylation, methylation_tmp])
    except Exception as e:
        print(f"Error reading index file: {e}")
        return None
            
    
    return methylation
    
def worker_function(task, motifs, counter, lock):
    """
    
    """
    contig, contig_length, mod_type, index_file = task
    
    try:
        result = calculate_contig_methylation_pattern(
            contig = contig, 
            contig_length=contig_length,
            motifs=motifs,
            mod_type=mod_type,
            index_file=index_file
        )
        with lock:
            counter.value += 1
        return result
    except:
        with lock:
            counter.value += 1
        return None


def data_split_methylation_parallel(contig_lengths, motifs, motif_index_dir, threads=1):
    """
    Calculate methylation pattern for each contig in the data split in parallel.
    """
    # Create and filter tasks: compile only those tasks with existing index files for each contig and modification type.
    tasks = [(contig, contig_lengths[contig], mod_type, os.path.join(motif_index_dir, f"{contig}_{mod_type}_motifs_positions.npz")) for contig in contig_lengths for mod_type in motifs if os.path.exists(os.path.join(motif_index_dir, f"{contig}_{mod_type}_motifs_positions.npz"))]

    # Create a progress manager
    manager = multiprocessing.Manager()
    counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=threads)

    # Create a process for the progress bar
    progress_bar_process = multiprocessing.Process(target=update_progress_bar, args=(counter, len(tasks), True))
    progress_bar_process.start()

    # Put them workers to work
    results = pool.starmap(worker_function, [(
        task, 
        motifs,
        counter,
        lock
        ) for task in tasks])
    results = [result for result in results if result is not None] #TODO: Check if this is necessary

    # Close the pool
    pool.close()
    pool.join()

    # Close the progress bar
    progress_bar_process.join()
    
    methylation_pattern = pl.concat(results)
    methylation_pattern = methylation_pattern\
            .sort(["mod_type", "motif", "contig"])
    
    return methylation_pattern
    

def create_methylation_matrix(methylation_features, motifs=None, min_motif_observations = 8):
    """
    Creates a feature matrix with methylation from motifs-scored or methylation features.
    """
    # check if the methylation features have the required columns
    required_columns = ["contig", "motif",  "mod_type", "mod_position", "n_mod", "n_nomod"]
    if not all(col in methylation_features.columns for col in required_columns):
        raise ValueError(f"Missing required columns in methylation features. Required columns: {', '.join(required_columns)}")
    
    # Calculate mean methylation for each motif
    matrix = methylation_features\
        .with_columns(
            mean = pl.col("n_mod") / (pl.col("n_mod") + pl.col("n_nomod")),
            motif_mod = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String),
            n_motifs = pl.col("n_mod") + pl.col("n_nomod")
        )\
        .filter(pl.col("n_motifs") >= min_motif_observations)
    
    if motifs:
        matrix = matrix.filter(pl.col("motif_mod").is_in(motifs))
    
    matrix = matrix.select(["contig", "motif_mod", "mean"])\
        .pivot(
            index = "contig",
            columns = "motif_mod",
            values = "mean",
            aggregate_function = "first"
        )\
        .fill_null(0.0)
    
    return matrix


def add_must_links(data, data_split, must_links):
    """
    Processes the must_links file and concatenates the filtered data to the main DataFrame.
    
    Parameters:
    data (pl.DataFrame): The main data DataFrame.
    must_links (pl.DataFrame): DataFrame containing must link pairs with columns ['ml_1', 'ml_2'].
    
    Returns:
    pl.DataFrame: Updated data DataFrame with concatenated must link pairs.
    """
    for ml_1, ml_2 in zip(must_links["ml_1"], must_links["ml_2"]):
        ml_1_data = data.filter(pl.col("") == ml_1).drop(pl.selectors.matches(".*(mean|var).*"))
        ml_2_data = data.filter(pl.col("") == ml_2).drop(pl.selectors.matches(".*(mean|var).*"))
        
        
        if ml_1_data.shape[0] == 0 or ml_2_data.shape[0] == 0:
            print(f"Must link {ml_1} or {ml_2} not found in data")
            continue
        
        assert ml_1_data.shape[0] == 1, f"Must link {ml_1} should have only one row."
        assert ml_2_data.shape[0] == 1, f"Must link {ml_2} should have only one row."
        ml = pl.concat([ml_1_data, ml_2_data], rechunk = True)
        data_split.extend(ml)
    
    return data_split


def generate_methylation_features(logger, args):
    
    
    paths = [args.motifs_scored, args.data, args.data_split, args.contig_fasta, args.bin_motifs]
    directories = [args.motif_index_dir]
    # if args.must_links:
    #     paths = paths + [args.must_links]
    check_files_exist(paths, directories)
    
    # check if output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Load the data
    motifs_scored, data, data_split, bin_consensus = load_data(args)
    
    # Get the unique contigs from the data split
    contigs = get_contigs(data_split)
    
    # Load the assembly file
    contig_sequences = read_fasta(args.contig_fasta, contigs)
    
    # Get lenths of the contigs
    contig_lengths = {contig: len(sequence) for contig, sequence in contig_sequences.items()}
    
    # Get the unique motifs
    motifs = get_motifs(
        motifs_scored = motifs_scored, 
        bin_consensus = bin_consensus,
        occurence_cutoff = args.motif_occurence_cutoff,
        min_motif_observations=args.min_motif_observations,
        ambiguous_interval = args.ambiguous_interval,
        ambiguous_motif_percentage_cutoff = args.ambiguous_motif_percentage_cutoff
    )
    
    if len(motifs) == 0:
        print(f"No motifs found with --motif-occurence-cutoff {args.motif_occurence_cutoff}, --min-motif-observations {args.min_motif_observations}, --ambiguous-interval {args.ambiguous_interval}, and --.")
        sys.exit(1)
    
    # Create methylation matrix for contig_splits
    print("Calculating methylation pattern for each contig split using multiple threads.")
    start_time = time.time()
    contig_split_methylation = data_split_methylation_parallel(contig_lengths, motifs, args.motif_index_dir, threads=args.num_process)
    end_time = time.time()  # Record the end time after function execution
    print(f"Methylation Pattern Execution time: {end_time - start_time} seconds")  # Print the execution time
    
    data_split_methylation_matrix = create_methylation_matrix(
        methylation_features = contig_split_methylation,
        min_motif_observations = args.min_motif_observations
    )
    
    # extract motfis from the data
    motifs_in_contig_split = contig_split_methylation\
        .with_columns(
            motif = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String)
        )\
        .get_column("motif")\
        .unique()\
        .to_list()
    
    data_split = data_split\
        .join(
            data_split_methylation_matrix,
            on = "contig",
            how = "left"
        )\
        .rename({"contig": ''})\
        .fill_nan(0.0)\
        .fill_null(0.0)
        
    # Create methylation matrix for all data
    motifs_scored_matrix = create_methylation_matrix(
        methylation_features = motifs_scored, 
        motifs = motifs_in_contig_split,
        min_motif_observations = args.min_motif_observations
    ).select(data_split_methylation_matrix.columns)
    
    data = data\
        .join(
            motifs_scored_matrix,
            on = "contig",
            how = "left"
        )\
        .rename({"contig": ''})\
        .fill_null(0.0)\
        .fill_nan(0.0)
    
    # If must_links are provided, add them to the data
    # if args.must_links:
    #     must_links = pl.read_csv(args.must_links, has_header = False, new_columns = ["ml_1", "ml_2"])
    #     data_split = add_must_links(data, data_split, must_links)
    
    
    
    try:
        data_split.write_csv(os.path.join(args.output, "data_split_methylation_matrix.csv"), separator=",", quote_style='never') 
        data.write_csv(os.path.join(args.output, "data_methylation_matrix.csv"), separator=",", quote_style='never')
    except Exception as e:
        print(f"An error occurred while writing the output: {e}")
        sys.exit(1)
    
    
    
