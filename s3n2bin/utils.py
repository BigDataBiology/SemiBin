import numpy as np
import os
import subprocess
from atomicwrites import atomic_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import pandas as pd
import random

def validate_args(args):

    def expect_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(
                    f"Error: Expected file '{f}' does not exist\n")
                sys.exit(1)

    def expect_file_list(fs):
        if fs is not None:
            for f in fs:
                expect_file(f)

    if args.cmd == 'predict_taxonomy':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)

    if args.cmd == 'generate_data_single':
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'generate_data_multi':
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'bin':
        expect_file(args.contig_fasta)
        expect_file(args.data)
        expect_file(args.data_split)

    if args.cmd == 'single_easy_bin':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

    if args.cmd == 'multi_easy_bin':
        if args.GTDB_reference is not None:
            expect_file(args.GTDB_reference)
        expect_file(args.contig_fasta)
        expect_file_list(args.bams)

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
    return np.clip(threshold, 4000, None)


def parse_mmseqs(mmseqs_result):
    species_result = mmseqs_result.query('rank_name == "species" and score > 0.95').values
    genus_result = mmseqs_result.query('rank_name == "genus" and score > 0.80').values

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


def generate_cannot_link(mmseqs_path,namelist,num_threshold,output,sample):
    import itertools
    mmseqs_result = pd.read_csv(mmseqs_path, sep='\t', header=None)
    mmseqs_result.columns = ['contig_name', 'taxon_ID', 'rank_name', 'scientific_name', 'temp_1', 'temp_2', 'temp_3',
                             'score', 'lineage']
    mmseqs_result = mmseqs_result[[
        'contig_name', 'rank_name', 'scientific_name', 'score', 'lineage']]
    mmseqs_result = mmseqs_result[mmseqs_result['contig_name'].isin(
        namelist)]
    cannot_link_species, cannot_link_genus, cannot_link_mix = parse_mmseqs(
        mmseqs_result)
    num_whole_data = np.clip(1000 * num_threshold, None, 4_000_000)
    max_num_mix = int(num_whole_data / 8)
    max_num_genus = int(max_num_mix / 2)
    max_num_species = num_whole_data - max_num_mix - max_num_genus
    if len(cannot_link_mix) > max_num_mix:
        cannot_link_mix = random.sample(cannot_link_mix, max_num_mix)

    if len(cannot_link_genus) > max_num_genus:
        cannot_link_genus = random.sample(cannot_link_genus, max_num_genus)

    if len(cannot_link_species) > max_num_species:
        cannot_link_species = random.sample(cannot_link_species, max_num_species)

    with open(os.path.join(output, f'{sample}.txt'), 'w') as out_cannot_link:
        for cannot in itertools.chain(
                cannot_link_species,
                cannot_link_genus,
                cannot_link_mix):
            out_cannot_link.write(f'{cannot[0]},{cannot[1]}\n')

def cal_num_bins(fasta_path, contig_output, hmm_output,
                 seed_output, binned_short,num_process):
    if not os.path.exists(contig_output + '.faa'):
        frag_out_log = open(contig_output + '.out', 'w')
        subprocess.check_call(
            ['run_FragGeneScan.pl',
             '-genome={}'.format(fasta_path),
             '-out={}'.format(contig_output),
             '-complete=0',
             '-train=complete',
             '-thread={}'.format(num_process),
             ],
            stdout=frag_out_log,
            stderr=subprocess.DEVNULL,
        )

    if not os.path.exists(hmm_output):
        hmm_out_log = open(hmm_output + '.out', 'w')
        subprocess.check_call(
            ['hmmsearch',
             '--domtblout',
             hmm_output,
             '--cut_tc',
             '--cpu', str(num_process),
             os.path.split(__file__)[0] + '/marker.hmm',
             contig_output + '.faa',
             ],
            stdout=hmm_out_log,
            stderr=subprocess.DEVNULL,
        )

    if not os.path.exists(seed_output):
        getmarker = os.path.split(__file__)[0] + '/test_getmarker.pl'
        subprocess.check_call(
            ['perl', getmarker,
             hmm_output,
             fasta_path,
             ('1001' if binned_short else '2501'), # threshold
             seed_output,
             ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def write_bins(namelist, contig_labels, output, contig_dict,
               recluster=False, origin_label=0,minfasta = 200000):
    from collections import defaultdict
    res = defaultdict(list)
    for label, name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    for label in res:
        bin = []
        whole_bin_bp = 0
        for contig in res[label]:
            rec = SeqRecord(
                Seq(str(contig_dict[contig])), id=contig, description='')
            bin.append(rec)
            whole_bin_bp += len(str(contig_dict[contig]))

        ofname = f'bin.{label}.fa' if not recluster \
                    else f'recluster_{origin_label}.bin.{label}.fa'
        if whole_bin_bp >= minfasta:
            with atomic_write(os.path.join(output, ofname), overwrite=True) as ofile:
                SeqIO.write(bin, ofile, 'fasta')


def cal_kl(m1, m2, v1, v2):
    m1 = np.clip(m1, 1e-6, None)
    m2 = np.clip(m2, 1e-6, None)
    v1 = np.clip(v1, 1.0, None)
    v2 = np.clip(v2, 1.0, None)
    value = np.log(np.sqrt(v2 / v1)) + \
        np.divide(np.add(v1, np.square(m1 - m2)), 2 * v2) - 0.5
    return np.clip(value, 1e-6, 1 - 1e-6)