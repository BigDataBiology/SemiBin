import numpy as np
import pandas as pd

def calculate_coverage(depth_file,threshold,edge=75,is_combined = False,contig_threshold=1000):
    """
    Input is position depth file generated from mosdepth or bedtools genomecov
    """
    data_contig_list = []
    coverage = []
    var = []

    split_contig_list = []
    split_coverage = []
    print('start')
    data = pd.read_csv(depth_file, sep='\t', header=None)
    data.columns = ['contig', 'start', 'end', 'value']
    print(data)
    groupby = data.groupby('contig')
    print(len(groupby))
    num = 0
    for contig, contig_data in groupby:
        num += 1
        if num % 50000 == 0:
            print(num)
        contig_data = contig_data.values
        depth_value = []
        if contig_data[len(contig_data) - 1][2] <= contig_threshold:
            continue
        if len(contig_data) == 1:
            if is_combined:
                coverage.append(contig_data[0][3])
                data_contig_list.append(contig)
                split_contig_list.append(contig + '_1')
                split_contig_list.append(contig + '_2')
                split_coverage.append(contig_data[0][3])
                split_coverage.append(contig_data[0][3])
                continue
            else:
                coverage.append(contig_data[0][3])
                var.append(0)
                data_contig_list.append(contig)
                continue

        for i in range(len(contig_data)):
            depth_value.extend([contig_data[i][3]] * (contig_data[i][2] - contig_data[i][1]))

        # not consider the edges of the contig
        depth_value_ = depth_value[edge:len(depth_value) - edge]
        mean = np.mean(depth_value_)
        coverage.append(mean)
        var.append(np.var(depth_value_))

        data_contig_list.append(contig)

        # generate coverage for must link pair
        if is_combined:
            if len(depth_value) >= threshold:
                depth_value1 = depth_value[0:len(depth_value) // 2]
                depth_value2 = depth_value[len(depth_value) // 2: len(depth_value)]
                split_contig_list.append(contig + '_1')
                split_contig_list.append(contig + '_2')
                depth_value1 = depth_value1[edge:len(depth_value1)]
                mean = np.mean(depth_value1)
                split_coverage.append(mean)
                depth_value2 = depth_value2[0:len(depth_value2) - edge]
                mean = np.mean(depth_value2)
                split_coverage.append(mean)

    if is_combined:
        data_contig_list = np.expand_dims(np.array(data_contig_list),axis=1)
        coverage = np.expand_dims(np.array(coverage),axis=1)
        contig_cov = np.concatenate((data_contig_list, coverage),axis=1)
        contig_cov = pd.DataFrame(data = contig_cov[:,1:],index = contig_cov[:,0],columns = ['cov'])
        contig_cov['cov'] = contig_cov['cov'].astype('float')
        split_contig_list = np.expand_dims(np.array(split_contig_list),axis=1)
        split_coverage = np.expand_dims(np.array(split_coverage),axis=1)
        split_contig_cov = np.concatenate((split_contig_list, split_coverage), axis=1)
        split_contig_cov = pd.DataFrame(data = split_contig_cov[:,1:],index = split_contig_cov[:,0],columns = ['cov'])
        split_contig_cov['cov'] = split_contig_cov['cov'].astype('float')
        return contig_cov , split_contig_cov
    else:
        data_contig_list = np.expand_dims(np.array(data_contig_list), axis=1)
        coverage = np.expand_dims(np.array(coverage), axis=1)
        var = np.expand_dims(np.array(var), axis=1)
        contig_cov = np.concatenate((data_contig_list, coverage), axis=1)
        contig_cov = np.concatenate((contig_cov, var), axis=1)
        contig_cov = pd.DataFrame(data=contig_cov[:, 1:], index=contig_cov[:, 0], columns=['mean', 'var'])
        contig_cov['mean'] = contig_cov['mean'].astype('float')
        contig_cov['var'] = contig_cov['var'].astype('float')
        return contig_cov
