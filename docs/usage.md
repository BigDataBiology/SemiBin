# Usage

## Generate cannot-link constrains

You can use  [mmseqs](https://github.com/soedinglab/MMseqs2) or [CAT](https://github.com/dutilh/CAT)(or other contig annotation tools) to get the taxonomic classifications of contigs. Then you can use the script `script/concatenate.py` to generate the cannot-link file(Format: contig_1,contig_2 ) that can be used in S<sup>3</sup>N<sup>2</sup>Bin. 

Contig annotation with mmseqs(GTDB reference genomes)

```bash
mmseqs createdb contig.fasta contig_DB

mmseqs taxonomy contig_DB GTDB taxonomyResult tmp

mmseqs createtsv contig_DB taxonomyResult taxonomyResult.tsv
```

Generate cannot-link constrains

```bash
python script/concatenate.py -i taxonomyResult.tsv -c contig.fasta -s sample-name -o output --mmseqs
```

Contig annotation with CAT

```bash
CAT contigs -c contig.fasta -d CAT_prepare_20200304/2020-03-04_CAT_database --path_to_prodigal path_to_prodigal --path_to_diamond path_to_diamond -t CAT_prepare_20200304/2020-03-04_taxonomy -o CAT_output/CAT --force -f 0.5 --top 11 --I_know_what_Im_doing --index_chunks 1

CAT add_names CAT_output/CAT.contig2classification.txt -o CAT_output/CAT.out -t CAT_prepare_20200304/2020-03-04_taxonomy --force --only_official
```

Generate cannot-link constrains

```bash
python script/concatenate.py -i CAT.out -c contig.fasta -s sample-name -o output --CAT
```

## Examples

### Single sample/co-assembly binning

#### Single sample/co-assembly binning pipeline

(1) Mapping reads to the contig fasta file. 

Bowtie2(Or any reads mapping tool)

```bash
bowtie2-build -f contig.fna contig.fna -p 16

bowtie2 -q --fr contig.fna -1 reads_1.fq.gz -2 reads_2.fq.gz -S contig.sam -p 64

samtools view -h -b -S contig.sam -o contig.bam -@ 64

samtools view -b -F 4 contig.bam -o contig.mapped.bam -@ 64

samtools sort -m 1000000000 contig.mapped.bam -o contig.mapped.sorted.bam -@ 64

samtools index contig.mapped.sorted.bam
```

(2) Generate cannot-link files for the contig fasta file.

(3) Run S<sup>3</sup>N<sup>2</sup>Bin.

```bash
S3N2Bin -i contig.fna -b *.bam -c cannot-link.txt -o output 
```

### Multi-samples binning(Must set -s parameter)

#### Multi-samples binning pipeline

(1) Concatenate all contigs from all samples together. Make sure that names of samples are unique and id for every contig is <sample_name><\separator><contig_id>. Concatenated contig format is:

```bash
<S1>O<Contig_1>
ATGCAAAA
<S1>O<Contig_2>
ATGCAAAA
<S1>O<Contig_3>
ATGCAAAA
<S2>O<Contig_1>
ATGCAAAA
<S2>O<Contig_2>
ATGCAAAA
<S3>O<Contig_1>
ATGCAAAA
<S3>O<Contig_2>
ATGCAAAA
<S3>O<Contig_3>
ATGCAAAA
```

(2) Map reads to the concatenated contig file to get the bam files.

(3) Generate cannot-link files for every sample. The name of the cannot-link file is <sample_name>.txt. Make sure the sample name here is same to that in step(1).

(4) Run S<sup>3</sup>N<sup>2</sup>Bin.

```bash
S3N2Bin -i whole_contig.fna -b *.bam -c *.txt -s C -o output
```

