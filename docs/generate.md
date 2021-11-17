# Generating the inputs to SemiBin

Starting with a metagenome, you need to generate a contigs file (`contigs.fna`)
and a sorted BAM file (`output.bam`).

1. Assemble it into a contigs FASTA file. In this case, we are using
   [NGLess](https://ngless.embl.de/) to combine FastQ preprocessing &amp;
   assembly (using
   [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884),
   but any other system will work.

```python
ngless "1.2"

input = paired('reads_1.fq.gz', 'reads_2.fq.gz')
input = preprocess(input) using |r|:
    r = substrim(r, min_quality=25)
    if len(r) < 45:
        discard

contigs = assemble(input)
write(contigs, ofile='contig.fna')
```

2. Map reads to the FASTA file

### Mapping using NGLess

```python
ngless "1.1"
import "samtools" version "1.0"

input = fastq('sample.fq.gz')
mapped = map(input, fafile='expected.fna')

write(samtools_sort(mapped),
    ofile='output.bam')
```

### Mapping using bowtie2

```bash
bowtie2-build -f contig.fna contig.fna -p 16

bowtie2 -q --fr -x contig.fna -1 reads_1.fq.gz -2 reads_2.fq.gz -S contig.sam -p 64

samtools view -h -b -S contig.sam -o contig.bam -@ 64

samtools view -b -F 4 contig.bam -o contig.mapped.bam -@ 64

samtools sort -m 1000000000 contig.mapped.bam -o contig.mapped.sorted.bam -@ 64

samtools index contig.mapped.sorted.bam
```

### Generate cannot-link constraints

You can also use [CAT](https://github.com/dutilh/CAT) for generating  contig taxonomic classifications and generating cannot-link file. See below for format
specifications if you want to use CAT.

```bash
CAT contigs \
        -c contig.fasta \
        -d CAT_prepare_20200304/2020-03-04_CAT_database \
        --path_to_prodigal path_to_prodigal \
        --path_to_diamond path_to_diamond \
        -t CAT_prepare_20200304/2020-03-04_taxonomy \
        -o CAT_output/CAT \
        --force \
        -f 0.5 \
        --top 11 \
        --I_know_what_Im_doing \
        --index_chunks 1

CAT add_names \
    CAT_output/CAT.contig2classification.txt \
    -o CAT_output/CAT.out \
    -t CAT_prepare_20200304/2020-03-04_taxonomy \
    --force \
    --only_official
```

Generate cannot-link constrains

```bash
python script/concatenate.py -i CAT.out -c contig.fasta -s sample-name -o output --CAT
```

