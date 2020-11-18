import numpy as np
import pandas as pd

def calculate_coverage(depth_file,threshold,edge = 75):
    """
    Input is position depth file generated from mosdepth or bedtools genomecov
    """
    data_contig_list = []
    coverage = []

    split_contig_list = []
    split_coverage = []

    data = pd.read_csv(depth_file, sep='\t', header=None)
    data.columns = ['contig', 'start', 'end', 'value']
    groupby = data.groupby('contig')

    for contig, contig_data in groupby:
        contig_data = contig_data.values
        depth_value = []
        for i in range(len(contig_data)):
            depth_value.extend([contig_data[i][3]] * (contig_data[i][2] - contig_data[i][1]))

        if len(depth_value) <= 1000:
            continue

        # not consider the edges of the contig
        depth_value_ = depth_value[edge:len(depth_value) - edge]
        mean = np.mean(depth_value_)
        if mean > 1:
            coverage.append(mean)
        else:
            coverage.append(0)
        data_contig_list.append(contig)

        # generate coverage for must link pair
        if len(depth_value) >= threshold:
            depth_value1 = depth_value[0:len(depth_value) // 2]
            depth_value2 = depth_value[len(depth_value) // 2: len(depth_value)]
            split_contig_list.append(contig + '_1')
            split_contig_list.append(contig + '_2')
            depth_value1 = depth_value1[edge:len(depth_value1)]
            mean = np.mean(depth_value1)
            if mean > 1:
                split_coverage.append(mean)
            else:
                split_coverage.append(0)

            depth_value2 = depth_value2[0:len(depth_value2) - edge]
            mean = np.mean(depth_value2)
            if mean > 1:
                split_coverage.append(mean)
            else:
                split_coverage.append(0)

    data_contig_list = np.array(data_contig_list).reshape(len(data_contig_list), 1)
    coverage = np.array(coverage).reshape(len(coverage), 1)
    split_contig_list = np.array(split_contig_list).reshape(len(split_contig_list), 1)
    split_coverage = np.array(split_coverage).reshape(len(split_coverage), 1)

    contig_cov = pd.DataFrame(np.concatenate((data_contig_list, coverage),axis=1))
    contig_cov = contig_cov.set_index(0)
    contig_cov.index.name = None

    split_contig_cov = pd.DataFrame(np.concatenate((split_contig_list, split_coverage), axis=1))
    split_contig_cov = split_contig_cov.set_index(0)
    split_contig_cov.index.name = None

    return contig_cov , split_contig_cov
