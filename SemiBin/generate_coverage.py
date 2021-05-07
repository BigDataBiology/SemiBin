import numpy as np
import pandas as pd


def calculate_coverage(depth_file, threshold, edge=75, is_combined=False,
                       contig_threshold=1000, sep=None, binned_thre_dict=None):
    """
    Input is position depth file generated from mosdepth or bedtools genomecov
    """
    data_contig_list = []
    coverage = []
    var = []

    split_contig_list = []
    split_coverage = []
    f = open(depth_file)

    def process(depth_value, contig_name):
        depth_value_ = depth_value[edge:len(depth_value) - edge]
        mean = np.mean(depth_value_)
        coverage.append(mean)
        var.append(np.var(depth_value_))
        data_contig_list.append(contig_name)
        if is_combined:
            if len(depth_value) >= threshold:
                depth_value1 = depth_value[0:len(depth_value) // 2]
                depth_value2 = depth_value[len(
                    depth_value) // 2: len(depth_value)]
                split_contig_list.append(contig_name + '_1')
                split_contig_list.append(contig_name + '_2')
                depth_value1 = depth_value1[edge:len(depth_value1)]
                mean = np.mean(depth_value1)
                split_coverage.append(mean)
                depth_value2 = depth_value2[0:len(depth_value2) - edge]
                mean = np.mean(depth_value2)
                split_coverage.append(mean)

    depth_value = []
    depth_contig = None

    for line in f:
        line_split = line.strip().split('\t')
        contig_name = line_split[0]
        length = int(line_split[2]) - int(line_split[1])
        value = int(line_split[3])

        if depth_contig is None:
            depth_contig = contig_name
        if contig_name == depth_contig:
            depth_value.extend([value] * length)
        if contig_name != depth_contig:
            if sep is None:
                cov_threshold = contig_threshold
            else:
                sample_name = depth_contig.split(sep)[0]
                cov_threshold = 1000 if sample_name in binned_thre_dict[1000] else 2500

            if len(depth_value) <= cov_threshold:
                depth_contig = contig_name
                depth_value = []
                depth_value.extend([value] * length)
                continue

            process(depth_value, depth_contig)
            depth_contig = contig_name
            depth_value = []
            depth_value.extend([value] * length)

    if depth_value != []:
        if sep is None:
            cov_threshold = contig_threshold
        else:
            sample_name = depth_contig.split(sep)[0]
            cov_threshold = 1000 if sample_name in binned_thre_dict[1000] else 2500
        if len(depth_value) > cov_threshold:
            process(depth_value, depth_contig)

    if is_combined:
        data_contig_list = np.expand_dims(np.array(data_contig_list), axis=1)
        coverage = np.expand_dims(np.array(coverage), axis=1)
        contig_cov = np.concatenate((data_contig_list, coverage), axis=1)
        contig_cov = pd.DataFrame(
            data=contig_cov[:, 1:], index=contig_cov[:, 0], columns=['cov'])
        contig_cov['cov'] = contig_cov['cov'].astype('float')
        split_contig_list = np.expand_dims(np.array(split_contig_list), axis=1)
        split_coverage = np.expand_dims(np.array(split_coverage), axis=1)
        split_contig_cov = np.concatenate(
            (split_contig_list, split_coverage), axis=1)
        split_contig_cov = pd.DataFrame(
            data=split_contig_cov[:, 1:], index=split_contig_cov[:, 0], columns=['cov'])
        split_contig_cov['cov'] = split_contig_cov['cov'].astype('float')
        return contig_cov, split_contig_cov
    else:
        data_contig_list = np.expand_dims(np.array(data_contig_list), axis=1)
        coverage = np.expand_dims(np.array(coverage), axis=1)
        var = np.expand_dims(np.array(var), axis=1)
        contig_cov = np.concatenate((data_contig_list, coverage), axis=1)
        contig_cov = np.concatenate((contig_cov, var), axis=1)
        contig_cov = pd.DataFrame(
            data=contig_cov[:, 1:], index=contig_cov[:, 0], columns=['mean', 'var'])
        contig_cov['mean'] = contig_cov['mean'].astype('float')
        contig_cov['var'] = contig_cov['var'].astype('float')
        return contig_cov
