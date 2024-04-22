import nanomotif as nm
from nanomotif.candidate import Motif
from collections import defaultdict
import re
import numpy as np
from multiprocessing import get_context
import os
import sys
import gzip
from Bio import SeqIO
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

def load_data(args, contigs):
    """Loads bed, assembly, and motifs files with robust error handling, exits on file not found."""
    try:
        # Check for the existence of the bed file
        if not os.path.exists(args.bed):
            raise FileNotFoundError(f"The bed file {args.bed} does not exist.")
        bed = nm.load_pileup(args.bed, min_fraction=args.min_fraction)

        # Check for the existence of the assembly file
        if not os.path.exists(args.assembly):
            raise FileNotFoundError(f"The assembly file {args.assembly} does not exist.")
        # assembly = nm.load_assembly(args.assembly)
        assembly = read_fasta(args.assembly, contigs)
        
        # Check for the existence of the motifs scored file
        if not os.path.exists(args.motifs_scored):
            raise FileNotFoundError(f"The motifs scored file {args.motifs_scored} does not exist.")
        motifs_scored = pl.read_csv(args.motifs_scored, separator="\t")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)  # Exit the program with an error code

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)  # Exit the program with an error code

    print("All files loaded successfully.")
    return bed, assembly, motifs_scored

# bed, assembly, motifs_scored = load_data(args)

def get_contigs(data_split):
    """
    Takes the data split and returns a list of unique contigs.
    """
    contigs_in_split = data_split["contig"].str.split("_").map_elements(lambda x: "_".join(x[:-1])).to_list()
    contigs_in_split = list(set(contigs_in_split))  # remove duplicates
    contigs_in_split.sort()  # Sort the list in place
    return contigs_in_split

def get_motifs(motifs_scored):
    """Extracts and returns unique motifs for each contig."""
    motifs = motifs_scored.select(["motif", "mod_position", "mod_type"]).unique()
    motifs_dict = motifs.to_dict()
    
    return [Motif(motif, position) for motif, position in zip(motifs_dict["motif"], motifs_dict["mod_position"])]


def subseq_indices(subseq, seq):
    """
    Find occurrence indices of a subseq in a seq, interpreting subseq using IUPAC codes
    """
    # Convert IUPAC subseq to regex pattern
    regex_pattern = ''.join(iupac_dict[char] for char in subseq)
    
    # Compile regex pattern
    compiled_subseq = re.compile(regex_pattern)
    
    # Find all matches using finditer
    index = [match.start() for match in re.finditer(compiled_subseq, seq)]
    
    return np.array(index)



def find_motif_indices(motifs, sequence):
    """Calculates forward and reverse indices for each motif in the given sequence."""
    motif_indices = defaultdict(lambda: {'forward': np.array([]), 'reverse': np.array([])})
    for motif in motifs:
        forward_indices = subseq_indices(motif.string, sequence) + motif.mod_position
        reverse_indices = subseq_indices(motif.reverse_compliment().string, sequence) + motif.mod_position
        # forward_indices = subseq_indices(motif.string, sequence.sequence) + motif.mod_position
        # reverse_indices = subseq_indices(motif.reverse_compliment().string, sequence.reverse_complement().sequence) + motif.mod_position
        motif_indices[f"{motif.string}_{motif.mod_position}"] = {'forward': forward_indices, 'reverse': reverse_indices}
    return motif_indices

def get_mod_type(motif):
    """
    Get the modification type from the motif.
    """
    seq, pos = motif.split("_")
    pos = int(pos)
    if seq[pos] == "C":
        mod_type = "m"
    elif seq[pos] == "A":
        mod_type = "a"
    
    return mod_type

def calculate_methylation_fraction(contig, bed, motif_indices, contig_length):
    """
    For each motif this will calculate the methylation fraction for each half of the contig.
    """
    contig_methylation = pl.DataFrame()
    
    for motif, indices in motif_indices.items():
        # Get mod_type from motif
        mod_type = get_mod_type(motif)
        
        pileup = bed.pileup.filter(
                pl.col("contig") == contig
            )\
            .filter(
                pl.col("mod_type") == mod_type
            )\
            .sort("position", descending = False)\
            .with_columns(
                pl.lit("").alias("motif")
            )
        
        
        # Get the actual NumPy arrays for forward and reverse indices
        forward_indices = pl.Series(indices['forward'])
        reverse_indices = pl.Series(indices['reverse'])

        # Apply motifs for the forward strand
        pileup = pileup.with_columns(
            pl.when(
                (pl.col("position").is_in(forward_indices)) & (pl.col("strand") == "+")
            ).then(pl.lit(motif)).otherwise(pl.col("motif")).alias("motif")
        )

        # Apply motifs for the reverse strand
        pileup = pileup.with_columns(
            pl.when(
                (pl.col("position").is_in(reverse_indices)) & (pl.col("strand") == "-")
            ).then(pl.lit(motif)).otherwise(pl.col("motif")).alias("motif")
        )
        
        
        # Divide pileup into half contig length
        pileup = pileup\
            .with_columns(
            pl.when(
                pl.col("position") < contig_length / 2
            ).then(pl.lit(contig + "_1")).otherwise(pl.lit(contig + "_2")).alias("contig")
        )
        ## Calculate the number of motifs in each half of the contig
        forward_indices_np = forward_indices.to_numpy()
        reverse_indices_np = reverse_indices.to_numpy()
        
        n_motifs_dict = {}
        n_motifs_dict[f"{contig}_1"] = len(forward_indices_np[forward_indices_np < contig_length / 2]) + len(reverse_indices_np[reverse_indices_np < contig_length / 2])
        n_motifs_dict[f"{contig}_2"] = len(forward_indices_np[forward_indices_np >= contig_length / 2]) + len(reverse_indices_np[reverse_indices_np >= contig_length / 2])
        
        motif_counts = pl.DataFrame({
            "motif": [motif, motif],
            "contig": [contig + "_1", contig + "_2"],
            "n_motifs": [n_motifs_dict[f"{contig}_1"], n_motifs_dict[f"{contig}_2"]]
        })
        
        # Filter out rows with empty motifs
        pileup = pileup.filter(pl.col("motif") != "")
        
        # Calculate fraction modified
        motif_mean_methylation = pileup\
            .with_columns(
                is_methylated = (pl.col("fraction_mod") >= 0.6).cast(pl.Int8)
            )\
            .group_by(["motif", "contig", "mod_type"])\
            .agg(
                pl.sum("is_methylated").alias("n_mod"),
            )\
            .join(
                motif_counts,
                on = ["motif", "contig"]
            )\
            .with_columns(
                n_nomod = pl.col("n_motifs") - pl.col("n_mod")
            )
        
        # Conditional check for cases where only one half of the contig has the motif or no motif at all
        if motif_mean_methylation.shape[0] < 2:
            if motif_mean_methylation.shape[0] == 1: # only one half has motif. Impute other half.
                
                contig_dict = {}
                contig_dict[f"{contig}_1"] = f"{contig}_2"
                contig_dict[f"{contig}_2"] = f"{contig}_1"
                
                contig_half_not_present = contig_dict[motif_mean_methylation["contig"][0]]
                
                motif_mean_methylation = pl.concat([motif_mean_methylation, pl.DataFrame({
                        "motif": [motif],
                        "contig": [contig_half_not_present],
                        "mod_type": [mod_type], # "m" or "a"
                        "n_mod": [None],
                        "n_motifs": [n_motifs_dict[contig_half_not_present]],
                        "n_nomod": [None]
                    })], rechunk = True)
            
            if motif_mean_methylation.shape[0] == 0:
                motif_mean_methylation = pl.DataFrame({
                    "motif": [motif, motif],
                    "contig": [contig + "_1", contig + "_2"],
                    "mod_type": [mod_type, mod_type],
                    "n_mod": [None, None],
                    "n_motifs": [n_motifs_dict[f"{contig}_1"], n_motifs_dict[f"{contig}_2"]],
                    "n_nomod": [None, None]
                }, schema = {
                    "motif": pl.String,
                    "contig": pl.String,
                    "mod_type": pl.String,
                    "n_mod": pl.Int64,
                    "n_motifs": pl.Int64,
                    "n_nomod": pl.Int64
                })
        
        contig_methylation = pl.concat([contig_methylation,motif_mean_methylation], rechunk = True)
        
    return contig_methylation




def worker_function(contig, bed, assembly, motifs):
    """
    Worker function that processes each contig in parallel.
    """
    # Get the sequence of the contig
    # contig_seq = assembly.assembly[contig]
    contig_seq = assembly[contig]
    
    # Find the motif indices in the contig sequence
    motif_indices = find_motif_indices(motifs, contig_seq)

    # Get the length of the contig
    # contig_length = len(contig_seq.sequence)
    contig_length = len(contig_seq)
    
    # Calculate the mean methylation for each motif in the contig halfs
    mean_methylation = calculate_methylation_fraction(contig, bed, motif_indices, contig_length)
    
    return mean_methylation


def process_contigs_parallel(data_split, args):
    # Get the unique contigs from the data split
    contigs = get_contigs(data_split)
    
    # Load the data
    bed, assembly, motifs_scored = load_data(args, contigs)
    
    # Check that all contigs are in the assembly and exit
    for contig in contigs:
        if contig not in assembly: #assembly.assembly:
            print(f"Error: Contig {contig} not found in the assembly.")
            sys.exit(1)
        
    # Get motifs in the motifs_scored
    motifs = get_motifs(motifs_scored)
    
    # Create a pool of workers
    pool = get_context("spawn").Pool(processes=getattr(args, 'threads', 1))

    # Put them workers to work
    results = pool.starmap(worker_function, [(
        contig, 
        bed, 
        assembly,
        motifs
        ) for contig in contigs])
    
    results = [result for result in results if result is not None]

    print(f"Processed {len(results)} contigs.")
    
    
    # Close the pool
    pool.close()
    pool.join()

    # Concatenate the results into a single DataFrame
    methylation_features = pl.concat(results, rechunk = True)
    
    methylation_features = methylation_features\
        .sort("motif", "contig", "mod_type")\
        .with_columns(
            [
                pl.col("motif")\
                .str.split_exact("_", 1)\
                .struct.rename_fields(["motif", "mod_position"])\
                .alias("fields"),
            ]
        )\
        .drop("motif")\
        .unnest("fields")
    
    return methylation_features

def create_methylation_matrix(methylation_features):
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
            motif_mod = pl.col("motif") + "_" + pl.col("mod_type") + "-" + pl.col("mod_position").cast(pl.String)
        )\
        .select(
            ["contig", "motif_mod", "mean"]
        )\
        .pivot(
            index = "contig",
            columns = "motif_mod",
            values = "mean",
            aggregate_function = "first"
        )\
        .fill_null("zero")
    
    # TODO: Impute missing values
    return matrix

