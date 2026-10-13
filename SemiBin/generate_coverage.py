from __future__ import annotations
import logging
import os
import subprocess
import shutil
import tempfile
from collections.abc import Iterable
from typing import TYPE_CHECKING
from .atomicwrite import atomic_write
from .utils import describe_returncode, format_stderr

if TYPE_CHECKING:
    import pandas as pd

def calculate_coverage(depth_stream: Iterable[str],
                       bam_file: str,
                       must_link_threshold: int,
                       edge: int = 75,
                       is_combined: bool = False,
                       contig_threshold: int = 1000,
                       sep: str | None = None,
                       sample_contig_threshold: dict[str, int] | None = None,
                       ) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Compute per-contig coverage statistics from per-base depth

    Parameters
    ----------
    depth_stream : iterable of str
        Lines in BED-graph format (contig, start, end, depth), such as the
        output of ``bedtools genomecov -bga``. Lines for the same contig must
        be consecutive.
    bam_file : str
        BAM filename (used to name the output columns and in error messages)
    must_link_threshold : int
        Contigs at least this long are split in two halves to generate
        must-link constraints (only used if ``is_combined``)
    edge : int, optional
        Number of bases to ignore at each end of the contig (coverage at the
        edges is less reliable)
    is_combined : bool, optional
        If True, compute only mean coverage (for both whole contigs and split
        contigs), to be combined with k-mer features. Otherwise, compute mean
        and variance of whole contigs.
    contig_threshold : int, optional
        Contigs shorter than this are ignored (used when ``sep`` is None)
    sep : str, optional
        Separator between sample name and contig name (multi-sample binning)
    sample_contig_threshold : dict[str, int], optional
        sample -> minimum contig length (used when ``sep`` is not None)

    Returns
    -------
    contig_cov : pd.DataFrame
        If ``is_combined``, one ``<bam_file>_cov`` column; otherwise
        ``<bam_file>_mean`` and ``<bam_file>_var`` columns
    split_contig_cov : pd.DataFrame or None
        Coverage of the split contigs (``<contig>_1`` and ``<contig>_2``) if
        ``is_combined``, otherwise None
    """
    import pandas as pd
    import numpy as np
    from itertools import groupby

    contigs = []
    mean_coverage = []
    var_coverage = []

    split_contigs = []
    split_coverage = []

    for contig_name, lines in groupby(depth_stream, lambda ell: ell.split('\t', 1)[0]):
        # Expand the run-length encoded BED-graph into per-base depth
        lengths = []
        depths = []
        for line in lines:
            line_split = line.strip().split('\t')
            length = int(line_split[2]) - int(line_split[1])
            depth = int(float(line_split[3]))
            lengths.append(length)
            depths.append(depth)
        depth_value = np.zeros(sum(lengths), dtype=int)
        pos = 0
        for length, depth in zip(lengths, depths, strict=True):
            depth_value[pos:pos+length] = depth
            pos += length

        if sep is None:
            cov_threshold = contig_threshold
        else:
            try:
                sample, contig = contig_name.split(sep)
                cov_threshold = sample_contig_threshold[sample]
            except ValueError:
                raise ValueError(f"Error parsing contig name '{contig_name}' in file {bam_file} (trying to split by {sep} separator)")
            except KeyError:
                raise KeyError(f"Error: {sample} not found (parsing {bam_file}, the original contig name was {contig_name} and using {sep} as separator)")
        if len(depth_value) < cov_threshold:
            continue
        inner_depth = depth_value[edge:-edge]
        mean_coverage.append(inner_depth.mean())
        var_coverage.append(inner_depth.var())
        contigs.append(contig_name)

        if is_combined:
            if len(depth_value) >= must_link_threshold:
                middle = len(inner_depth) // 2
                split_contigs.append(contig_name + '_1')
                split_contigs.append(contig_name + '_2')

                split_coverage.append(inner_depth[:middle].mean())
                split_coverage.append(inner_depth[middle:].mean())

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
            '{0}_var'.format(bam_file): var_coverage,
            }, index=contigs), None


def generate_cov(bam_file: str,
                 bam_index: int,
                 out: str,
                 must_link_threshold: int,
                 is_combined: bool,
                 contig_threshold: int | dict[str, int],
                 logger: logging.Logger,
                 sep: str | None = None) -> str:
    """
    Call bedtools and generate coverage file(s)

    Writes ``<bam basename>_<bam_index>_data_cov.csv`` (and, if
    ``is_combined``, ``<bam basename>_<bam_index>_data_split_cov.csv``) to
    ``out``.

    Parameters
    ----------
    bam_file : str
        Sorted BAM file
    bam_index : int
        Index of the BAM file in the input list (used to disambiguate output
        filenames if two BAM files have the same basename)
    out : str
        Output directory
    must_link_threshold : int
        Contigs at least this long are split to generate must-link constraints
    is_combined : bool
        Whether coverage is to be combined with k-mer features as input to the
        neural network (used when there are 5 or more samples)
    contig_threshold : int or dict[str, int]
        Minimum length of contigs to be binned. In multi-sample mode (``sep``
        is not None), a dictionary sample -> minimum length
    logger : logging.Logger
    sep : str, optional
        Separator for multi-sample binning

    Returns
    -------
    bam_file : str
        The input ``bam_file`` (so that callers can report progress)
    """
    import numpy as np
    logger.debug(f'Processing `{bam_file}`')
    bam_name = os.path.split(bam_file)[-1] + f'_{bam_index}'

    bedtools_path = shutil.which('bedtools')
    if bedtools_path is None:
        raise OSError('bedtools binary not found in PATH (it is required to compute coverage from BAM files). '
                      'Please install bedtools.')

    with tempfile.TemporaryFile() as bed_err:
        bed_p = subprocess.Popen(
            [bedtools_path, 'genomecov',
             '-bga',
             '-ibam', bam_file],
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=bed_err)

        try:
            contig_cov, must_link_contig_cov = calculate_coverage(
                                            bed_p.stdout,
                                            bam_file,
                                            must_link_threshold,
                                            is_combined=is_combined,
                                            sep=sep,
                                            contig_threshold=(contig_threshold if sep is None else 1000),
                                            sample_contig_threshold=(contig_threshold if sep is not None else None))
        except BaseException:
            bed_p.kill()
            raise
        finally:
            bed_p.stdout.close()
            rc = bed_p.wait()

        if rc != 0:
            bed_err.seek(0)
            err = format_stderr(bed_err.read())
            raise OSError(f"Running `bedtools genomecov` {describe_returncode(rc)} (processing {bam_file}).{err}")

    if len(contig_cov) == 0:
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
        with atomic_write(os.path.join(out, f'{bam_name}_data_cov.csv'), overwrite=True) as ofile:
            contig_cov.to_csv(ofile)

        must_link_contig_cov = must_link_contig_cov.apply(lambda x: x + 1e-5)
        if sep is None:
            abun_split_scale = (must_link_contig_cov.mean() / 100).apply(np.ceil) * 100
            must_link_contig_cov = must_link_contig_cov.div(abun_split_scale)

        with atomic_write(os.path.join(out, f'{bam_name}_data_split_cov.csv'), overwrite=True) as ofile:
            must_link_contig_cov.to_csv(ofile)

    return bam_file

def combine_cov(cov_dir: str,
                bam_list: list[str],
                is_combined: bool) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """
    Combine the per-BAM coverage files (written by ``generate_cov``) into a
    single table

    Parameters
    ----------
    cov_dir : str
        Directory where coverage files are stored
    bam_list : list of str
        BAM files (in the same order as passed to ``generate_cov``)
    is_combined : bool
        Whether to also process split files

    Returns
    -------
    cov : pd.DataFrame
        One column per BAM file
    split_cov : pd.DataFrame or None
        Coverage of split contigs if ``is_combined``, otherwise None
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

def generate_cov_from_abundances(
        abundances: list[str],
        output: str,
        contig_path: str,
        contig_threshold: int = 1000,
        sep: str | None = None,
        sample_contig_threshold: dict[str, int] | None = None,
        ) -> tuple[pd.DataFrame, pd.DataFrame] | pd.DataFrame:
    """Build coverage features from pre-computed abundance files

    The abundance files (e.g., from ``strobealign --aemb``) are expected to
    contain abundances of the split contigs (``<contig>_1`` and
    ``<contig>_2``); the abundance of the whole contig is computed as the
    average of the two halves.

    Parameters
    ----------
    abundances : list of str
        Abundance files (tab-separated, no header: contig name, abundance), one
        per sample
    output : str
        Output directory
    contig_path : str
        Contig FASTA file
    contig_threshold : int, optional
        Minimum length of contigs to be binned (used when ``sep`` is None)
    sep : str, optional
        Separator for multi-sample binning
    sample_contig_threshold : dict[str, int], optional
        sample -> minimum contig length (used when ``sep`` is not None)

    Returns
    -------
    If ``sep`` is None, a tuple ``(cov, cov_split)`` (which are also written
    to ``cov.csv`` and ``cov_split.csv`` in ``output``). Otherwise, only
    ``cov_split`` (unnormalized), which is to be split by sample by the caller.
    """
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
                try:
                    sample_name, contig_name = h.split(sep)
                    cov_threshold = sample_contig_threshold[sample_name]
                except ValueError:
                    raise ValueError(f"Error parsing contig name '{h}' in file {abun_file} (trying to split by {sep} separator)")
                except KeyError:
                    raise KeyError(f"Error: {sample_name} not found")

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
