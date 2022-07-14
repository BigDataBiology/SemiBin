import os
import subprocess
from atomicwrites import atomic_write

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
            }, index=contigs)


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
    logger.info('Processing `{}`'.format(bam_file))
    bam_name = os.path.split(bam_file)[-1] + '_{}'.format(bam_index)

    bed_p = subprocess.Popen(
        ['bedtools', 'genomecov',
         '-bga',
         '-ibam', bam_file],
        universal_newlines=True,
        stdout=subprocess.PIPE)

    if is_combined:
        contig_cov, must_link_contig_cov = calculate_coverage(bed_p.stdout, bam_file, threshold, is_combined = is_combined, sep = sep, contig_threshold = contig_threshold if sep is None else 1000, contig_threshold_dict =  contig_threshold if sep is not None else None)
        if bed_p.wait() != 0:
            raise OSError("Failure in running bedtools")

        contig_cov = contig_cov.apply(lambda x: x + 1e-5)
        must_link_contig_cov = must_link_contig_cov.apply(lambda x: x + 1e-5)
        if sep is None:
            abun_scale = (contig_cov.mean() / 100).apply(np.ceil) * 100
            abun_split_scale = (must_link_contig_cov.mean() / 100).apply(np.ceil) * 100
            contig_cov = contig_cov.div(abun_scale)
            must_link_contig_cov = must_link_contig_cov.div(abun_split_scale)

        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)

        with atomic_write(os.path.join(out, '{}_data_split_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            must_link_contig_cov.to_csv(ofile)
    else:
        contig_cov = calculate_coverage(bed_p.stdout, bam_file, threshold, is_combined=is_combined, sep = sep,
                                        contig_threshold = contig_threshold if sep is None else 1000,
                                        contig_threshold_dict =  contig_threshold if sep is not None else None)
        if bed_p.wait() != 0:
            raise OSError("Failure in running bedtools")

        contig_cov = contig_cov.apply(lambda x: x + 1e-5)

        with atomic_write(os.path.join(out, '{}_data_cov.csv'.format(bam_name)), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)

    return (bam_file, logger)

def combine_cov(cov_dir, bam_list, is_combined):
    """
    generate cov/cov_split for every sample in one file
    """
    import pandas as pd
    data_cov = pd.read_csv(os.path.join(cov_dir, '{}_data_cov.csv'.format(
        os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)
    if is_combined:
        data_split_cov = pd.read_csv(os.path.join(cov_dir, '{}_data_split_cov.csv'
                                    .format(os.path.split(bam_list[0])[-1] + '_{}'.format(0))), index_col=0)

    for bam_index, bam_file in enumerate(bam_list):
        if bam_index == 0:
            continue
        cov = pd.read_csv(os.path.join(cov_dir, '{}_data_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index))),index_col=0)
        cov.index = cov.index.astype(str)
        data_cov = pd.merge(data_cov, cov, how='inner', on=None,
                            left_index=True, right_index=True, sort=False, copy=True)

        if is_combined:
            cov_split = pd.read_csv(os.path.join(cov_dir, '{}_data_split_cov.csv'.format(
                os.path.split(bam_file)[-1] + '_{}'.format(bam_index))), index_col=0)

            data_split_cov = pd.merge(data_split_cov, cov_split, how='inner', on=None,
                                      left_index=True, right_index=True, sort=False, copy=True)
    data_cov.index = data_cov.index.astype(str)
    if is_combined:
        return data_cov, data_split_cov
    else:
        return data_cov

