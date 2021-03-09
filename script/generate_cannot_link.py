"""
A script to generate cannot-link constrains from the output of CAT.
"""
import argparse
import pandas as pd
from Bio import SeqIO
import math
import numpy as np
import random
import os


def generate_CAT(cat_result):
    must_link = []
    cannot_link_genus = []
    cannot_link_species = []

    num_genus = 0
    num_species = 0

    for i in range(len(cat_result)):
        if cat_result[i][1] == 'unclassified':
            continue
        genus1 = cat_result[i][2]

        if isinstance(genus1, float):
            if math.isnan(genus1):
                continue

        for j in range(i + 1, len(cat_result)):
            if cat_result[j][1] == 'unclassified':
                continue

            genus2 = cat_result[j][2]
            if isinstance(genus2, float):
                if math.isnan(genus2):
                    continue

            if genus1 == 'not classified' or genus2 == 'not classified':
                continue

            genus_score1 = genus1.split(':')[-1]
            genus_score2 = genus2.split(':')[-1]
            if float(genus_score1) <= 0.8 or float(genus_score2) <= 0.8:
                continue
            genus_name1 = genus1.strip(':' + genus_score1)
            genus_name2 = genus2.strip(':' + genus_score2)
            if genus_name1 != genus_name2:
                cannot_link_genus.append((cat_result[i][0], cat_result[j][0]))
                num_genus += 1
            else:
                species1 = cat_result[i][3]
                species2 = cat_result[j][3]
                if isinstance(species1, float) or isinstance(species2, float):
                    if math.isnan(species1) or math.isnan(species2):
                        continue

                if species1 == 'not classified' or species2 == 'not classified':
                    continue
                species_score1 = species1.split(':')[-1]
                species_score2 = species2.split(':')[-1]
                if float(species_score1) <= 0.95 or float(species_score2) <= 0.95:
                    continue
                species_name1 = species1.strip(':' + species_score1)
                species_name2 = species2.strip(':' + species_score2)
                if species_name1 != species_name2:
                    cannot_link_species.append(
                        (cat_result[i][0], cat_result[j][0]))
                    num_species += 1
                else:
                    must_link.append((cat_result[i][0], cat_result[j][0]))
    return must_link, cannot_link_genus, cannot_link_species


def generate_mmseqs(mmseqs_file):
    species_result = mmseqs_file[(mmseqs_file['rank_name'] == 'species') & (
        mmseqs_file['score'] > 0.95)].values
    genus_result = mmseqs_file[(mmseqs_file['rank_name'] == 'genus') & (
        mmseqs_file['score'] > 0.80)].values

    cannot_link_species = []
    for i in range(len(species_result)):
        for j in range(i + 1, len(species_result)):
            if species_result[i][2] != species_result[j][2]:
                cannot_link_species.append(
                    (species_result[i][0], species_result[j][0]))

    cannot_link_genus = []
    for i in range(len(genus_result)):
        for j in range(i + 1, len(genus_result)):
            if genus_result[i][2] != genus_result[j][2]:
                cannot_link_genus.append(
                    (genus_result[i][0], genus_result[j][0]))

    cannot_link_mix = []
    for i in range(len(species_result)):
        genus_name = species_result[i][4].split(';')[-2]
        for j in range(len(genus_result)):
            if genus_name != genus_result[j][2]:
                cannot_link_mix.append(
                    (species_result[i][0], genus_result[j][0]))

    return cannot_link_species, cannot_link_genus, cannot_link_mix


def get_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    basepair_sum = 0
    threshold = 0
    whole_len = np.sum(contig_len)
    contig_len.sort(reverse=True)
    index = 0
    while(basepair_sum / whole_len < 0.98):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    threshold = max(threshold, 4000)
    return threshold


def generate_file(annotation_file, contig_file, output, sample, tool=None):
    whole_contig_bp = 0
    contig_bp_2500 = 0
    contig_length_list = []
    for seq_record in SeqIO.parse(contig_file, "fasta"):
        if len(seq_record) >= 1000 and len(seq_record) <= 2500:
            contig_bp_2500 += len(seq_record)
        whole_contig_bp += len(seq_record)
        contig_length_list.append(len(seq_record))

    must_link_threshold = get_threshold(contig_length_list)
    binned_short = contig_bp_2500 / whole_contig_bp < 0.05
    threshold = 1000 if binned_short else 2500
    namelist = []
    num_threshold = 0

    for seq_record in SeqIO.parse(contig_file, "fasta"):
        if len(seq_record) > threshold:
            namelist.append(seq_record.id)
        if len(seq_record) >= must_link_threshold:
            num_threshold += 1
    os.makedirs(output, exist_ok=True)

    if tool == 'CAT':
        cat_result = pd.read_csv(annotation_file, sep='\t')
        cat_result = cat_result[['# contig',
                                 'classification', 'genus', 'species']]
        cat_result = cat_result[cat_result['# contig'].isin(namelist)].values

        _, cannot_link_genus, cannot_link_species = generate_CAT(
            cat_result)

        if len(cannot_link_genus) > 4000000:
            if num_threshold * 1000 < 4000000:
                cannot_link_genus = random.sample(
                    cannot_link_genus, num_threshold * 1000)
            else:
                cannot_link_genus = random.sample(cannot_link_genus, 4000000)
        else:
            if num_threshold * 1000 < len(cannot_link_genus):
                cannot_link_genus = random.sample(
                    cannot_link_genus, num_threshold * 1000)
        out_text = open(output + '/{}.txt'.format(sample), 'w')
        for cannot in cannot_link_species:
            out_text.write(cannot[0] + ',' + cannot[1])
            out_text.write('\n')

        for cannot in cannot_link_genus:
            out_text.write(cannot[0] + ',' + cannot[1])
            out_text.write('\n')

    if tool == 'mmseqs':
        mmseqs_result = pd.read_csv(annotation_file, sep='\t', header=None)
        mmseqs_result.columns = ['contig_name', 'taxon_ID', 'rank_name', 'scientific_name', 'temp_1', 'temp_2', 'temp_3',
                                 'score', 'lineage']
        mmseqs_result = mmseqs_result[[
            'contig_name', 'rank_name', 'scientific_name', 'score', 'lineage']]
        mmseqs_result = mmseqs_result[mmseqs_result['contig_name'].isin(
            namelist)]
        cannot_link_species, cannot_link_genus, cannot_link_mix = generate_mmseqs(
            mmseqs_result)
        num_whole_data = 1000 * num_threshold if 1000 * \
            num_threshold < 4000000 else 4000000
        num_mix = int(num_whole_data / 8)
        num_genus = int(num_mix / 2)
        num_species = num_whole_data - num_mix - num_genus
        if len(cannot_link_mix) > num_mix:
            cannot_link_mix = random.sample(cannot_link_mix, num_mix)

        if len(cannot_link_genus) > num_genus:
            cannot_link_genus = random.sample(cannot_link_genus, num_genus)

        if len(cannot_link_species) > num_species:
            cannot_link_species = random.sample(cannot_link_species, num_species)

        out_text = open(output + '/{}.txt'.format(sample), 'w')
        for cannot in cannot_link_species:
            out_text.write(cannot[0] + ',' + cannot[1])
            out_text.write('\n')

        for cannot in cannot_link_genus:
            out_text.write(cannot[0] + ',' + cannot[1])
            out_text.write('\n')

        for cannot in cannot_link_mix:
            out_text.write(cannot[0] + ',' + cannot[1])
            out_text.write('\n')




def main():
    parser = argparse.ArgumentParser(
        description="Generate cannt-link constrains from the output of CAT")
    parser.add_argument('-i', '--input-file',
                        required=True,
                        help='Path to the input CAT output files.',
                        dest='input_files',
                        default=None)
    parser.add_argument('-c', '--contig-file',
                        required=True,
                        help='Path to the contig fasta file corresponding to the CAT output file.',
                        dest='contig_file')
    parser.add_argument('-s', '--sample-name',
                        required=True,
                        help='Sample name used in the output file when multiple samples binning',
                        dest='sample',
                        default=None)
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None
                        )
    parser.add_argument('--CAT',
                        help='Input contig annotation from CAT',
                        action='store_true',
                        dest='CAT')
    parser.add_argument('--mmseqs',
                        help='Input contig annotation from mmseqs',
                        action='store_true',
                        dest='mmseqs')

    args = parser.parse_args()
    if args.CAT:
        generate_file(
            args.input_files,
            args.contig_file,
            args.output,
            args.sample,
            tool='CAT')
    if args.mmseqs:
        generate_file(
            args.input_files,
            args.contig_file,
            args.output,
            args.sample,
            tool='mmseqs')


if __name__ == '__main__':
    main()
