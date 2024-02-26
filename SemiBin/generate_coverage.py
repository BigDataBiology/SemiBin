import os
import subprocess
from .atomicwrite import atomic_write

def calculate_coverage(depth_stream, bam_file, must_link_threshold, edge=75, is_combined=False,
                       contig_threshold=1000, sep=None, contig_threshold_dict=None):
    """
    depth_stream : an iterable like the output of bedtools genomecov
    bam_file : str filename
    must_link_threshold : int
    edge : int
    is_combined : bool

    """
    import pandas as pd
    import numpy as np
    from itertools import groupby

    contigs = []
    mean_coverage = []
    var = []

    split_contigs = []
    split_coverage = []

    for contig_name, lines in groupby(depth_stream, lambda ell: ell.split('\t', 1)[0]):
        lengths = []
        values = []
        for line in lines:
            line_split = line.strip().split('\t')
            length = int(line_split[2]) - int(line_split[1])
            value = int(float(line_split[3]))
            lengths.append(length)
            values.append(value)
        depth_value = np.zeros(sum(lengths), dtype=int)
        s = 0
        for ell,v in zip(lengths, values):
            depth_value[s:s+ell] = v
            s += ell

        if sep is None:
            cov_threshold = contig_threshold
        else:
            sample_name = contig_name.split(sep)[0]
            cov_threshold = contig_threshold_dict[sample_name]
        if len(depth_value) < cov_threshold:
            continue
        depth_value_ = depth_value[edge:-edge]
        mean_coverage.append(depth_value_.mean())
        var.append(depth_value_.var())
        contigs.append(contig_name)

        if is_combined:
            if len(depth_value) >= must_link_threshold:
                middle = len(depth_value_) // 2
                split_contigs.append(contig_name + '_1')
                split_contigs.append(contig_name + '_2')

                split_coverage.append(depth_value_[:middle].mean())
                split_coverage.append(depth_value_[middle:].mean())

    if is_combined:
        contig_cov = pd.DataFrame(
                { '{0}_cov'.format(bam_file): mean_coverage,
                }, index=contigs)
        split_contig_cov = pd.DataFrame(
                { '{0}_cov'.format(bam_file): split_coverage,
                }, index=split_contigs)
        return contig_cov, split_contig_cov
    else:
        return pd.DataFrame({
            '{0}_mean'.format(bam_file): mean_coverage,
            '{0}_var'.format(bam_file): var,
            }, index=contigs), None


def generate_cov(bam_file, bam_index, out, threshold,
                 is_combined, contig_threshold, logger, sep = None):
    """
    Call bedtools and generate coverage file

    bam_file: bam files used
    out: output
    threshold: threshold of contigs that will be binned
    is_combined: if using abundance feature in deep learning. True: use
    contig_threshold: threshold of contigs for must-link constraints
    sep: separator for multi-sample binning
    """
    import numpy as np
    logger.debug('Processing `{}`'.format(bam_file))
    bam_name = os.path.split(bam_file)[-1] + '_{}'.format(bam_index)

    bed_p = subprocess.Popen(
        ['bedtools', 'genomecov',
         '-bga',
         '-ibam', bam_file],
        universal_newlines=True,
        stdout=subprocess.PIPE)

    contig_cov, must_link_contig_cov = calculate_coverage(
                                    bed_p.stdout,
                                    bam_file,
                                    threshold,
                                    is_combined=is_combined,
                                    sep=sep,
                                    contig_threshold=(contig_threshold if sep is None else 1000),
                                    contig_threshold_dict=(contig_threshold if sep is not None else None))

    if bed_p.wait() != 0:
        raise OSError(f"Failure in running bedtools ({bam_file})")
    elif len(contig_cov) == 0:
        logger.critical(f"Running `bedtools genomecov` did not return an error, but the result is an empty file (processing {bam_file}). "
                        "Please check your input files: SemiBin expects that they are sorted BAM files.")
        raise OSError(f"Running bedtools returned an empty file ({bam_file})")

    contig_cov += 1e-5
    if not is_combined:
        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)
    else:
        if sep is None:
            abun_scale = (contig_cov.mean() / 100).apply(np.ceil) * 100
            contig_cov = contig_cov.div(abun_scale)
        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)

    if is_combined:
        must_link_contig_cov = must_link_contig_cov.apply(lambda x: x + 1e-5)
        if sep is None:
            abun_split_scale = (must_link_contig_cov.mean() / 100).apply(np.ceil) * 100
            must_link_contig_cov = must_link_contig_cov.div(abun_split_scale)

        with atomic_write(os.path.join(out, '{}_data_split_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            must_link_contig_cov.to_csv(ofile)

    return bam_file

def combine_cov(cov_dir : str, bam_list, is_combined : bool): # bam_list : list[str] does not work on Python 3.7
    """
    generate cov/cov_split for every sample in one file

    Parameters
    ----------
    cov_dir : where coverage files are stored
    bam_list : list of BAM files
    is_combined : whether to process split files

    Returns
    -------
    cov : DataFrame
    split_cov : DataFrame (if is_combined) or None (otherwise)
    """
    import pandas as pd

    covs = []
    split_covs = []
    for bam_index, bam_file in enumerate(bam_list):
        bam_fname = os.path.split(bam_file)[-1]
        covs.append(
                pd.read_csv(f'{cov_dir}/{bam_fname}_{bam_index}_data_cov.csv', index_col=0))

        if is_combined:
            split_covs.append(
                    pd.read_csv(f'{cov_dir}/{bam_fname}_{bam_index}_data_split_cov.csv', index_col=0))

    data_cov = pd.concat(covs, axis=1)
    data_cov.index = data_cov.index.astype(str)
    if is_combined:
        data_split_cov = pd.concat(split_covs, axis=1)
        data_split_cov.index = data_split_cov.index.astype(str)
    else:
        data_split_cov = None
    return data_cov, data_split_cov

def generate_cov_from_abundances(abundances, output, contig_path, contig_threshold=1000, sep=None, contig_threshold_dict=None):
    import pandas as pd
    import numpy as np
    from .fasta import fasta_iter

    abun_split_list = []
    for abun_file in abundances:
        abun_split_data = pd.read_csv(abun_file, sep='\t', index_col=0, header=None)
        abun_split_data.index.name = None
        abun_split_data.columns = [f'{abun_file}']

        binned_contig = []
        for h, seq in fasta_iter(contig_path):
            if sep is None:
                cov_threshold = contig_threshold
            else:
                sample_name = h.split(sep)[0]
                cov_threshold = contig_threshold_dict[sample_name]

            if len(seq) >= cov_threshold:
                binned_contig.append(h + '_1')
                binned_contig.append(h + '_2')
        abun_split_data = abun_split_data.loc[binned_contig]
        abun_split_data = abun_split_data.apply(lambda x: x + 1e-5)
        if sep is None:
            abun_scale = (abun_split_data.mean() / 100).apply(np.ceil) * 100
            abun_split_data = abun_split_data.div(abun_scale)
        abun_split_list.append(abun_split_data)

    abun_split = pd.concat(abun_split_list, axis=1)
    if sep is None:
        with atomic_write(os.path.join(output, 'cov_split.csv'), overwrite=True) as ofile:
            abun_split.to_csv(ofile)

        index_name = abun_split.index.tolist()
        data_index_name = []
        for i in range(len(index_name) // 2):
            data_index_name.append(index_name[2 * i][0:-2])

        abun_split_values = abun_split.values
        abun = (abun_split_values[::2] + abun_split_values[1::2]) / 2
        columns = [f'{abun_file}' for abun_file in abundances]
        abun = pd.DataFrame(abun, index=data_index_name, columns=columns)
        with atomic_write(os.path.join(output, 'cov.csv'), overwrite=True) as ofile:
            abun.to_csv(ofile)
        return abun, abun_split
    else:
        return abun_split