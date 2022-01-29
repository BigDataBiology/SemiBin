# Usage

## Single-sample binning

Input: S1.fna and S1.bam

### Easy single binning mode

```bash
SemiBin single_easy_bin -i S1.fna -b S1.bam -o output 
```
Or with our built-in model(`human_gut`/`dog_gut`/`ocean`/`soil`/`cat_gut`/`human_oral`/`mouse_gut`/`pig_gut`/`built_environment`/`wastewater`/`whole`)
```bash
SemiBin single_easy_bin -i S1.fna -b S1.bam -o output --environment human_gut
```
### Advanced workflows

(1)  Generate `data.csv/data_split.csv` 
```bash
SemiBin generate_data_single -i S1.fna -b S1.bam -o S1_output
```
(2) Generate cannot-link 
```bash
SemiBin generate_cannot_links -i S1.fna -o S1_output
```
(3) Train
```bash
SemiBin train -i S1.fna --data S1_output/train.csv --data-split S1_output/train_split.csv -c S1_output/cannot/cannot.txt -o S1_output --mode single
```
(4) Bin 
```bash
SemiBin bin -i S1.fna --model S1_output/model.h5 --data S1_output/data.csv -o output
```
or with our built-in model(`human_gut`/`dog_gut`/`ocean`/`soil`/`cat_gut`/`human_oral`/`mouse_gut`/`pig_gut`/`built_environment`/`wastewater`/`whole`)
```bash
SemiBin bin -i S1.fna --data S1_output/data.csv -o output --environment human_gut
```

### SemiBin(pretrain)
You can train a model from several samples. 

If you have S1.fna, S1/data.csv,  S1/data_split.csv, S1/cannot/cannot.txt ; S2.fna, S2/data.csv,  S2/data_split.csv, S2/cannot/cannot.txt; S3.fna, S3/data.csv,  S3/data_split.csv, S3/cannot/cannot.txt. You can train the model from 3 samples.

```bash
SemiBin train -i S1.fna S2.fna S3.fna --data S1/train.csv S2/train.csv S3/train.csv --data-split S1/train_split.csv S2/train_split.csv S3/train_split.csv -c S1/cannot.txt s2/cannot.txt S3/cannot.txt -o output --mode several 
```

## Co-assembly binning

Input: contig.fna and S1.bam, S2.bam, S3.bam

### Easy single binning mode
```bash
SemiBin single_easy_bin -i contig.fna -b S1.bam S2.bam S3.bam -o output
```
### Advanced workflows

(1)  Generate `data.csv/data_split.csv` 
```bash
SemiBin generate_data_single -i contig.fna -b S1.bam S2.bam S3.bam -o contig_output
```
(2) Generate cannot-link 
```bash
SemiBin generate_cannot_links -i contig.fna -o contig_output
```
(3) Train
```bash
SemiBin train -i contig.fna --data contig_output/train.csv --data-split contig_output/train_split.csv -c contig_output/cannot/cannot.txt -o contig_output --mode single
```
(4) Bin
```bash
SemiBin bin -i contig.fna --model contig_output/model.h5 --data contig_output/data.csv -o output
```


## Multi-sample binning

Input: 
original fasta: S1.fna S2.fna S3.fna S4.fna S5.fna 
combined: combined.fna and S1.bam, S2.bam, S3.bam, S4.bam, S5.bam

The format of combined.fna: for every contig, format of the name is `<sample_name>:<contig_name>`, where
`:` is the default separator (it can be changed with the `--separator`
argument). *Note:* Make sure the sample names are unique and  the separator
does not introduce confusion when splitting. For example:

```bash
>S1:Contig_1
AGATAATAAAGATAATAATA
>S1:Contig_2
CGAATTTATCTCAAGAACAAGAAAA
>S1:Contig_3
AAAAAGAGAAAATTCAGAATTAGCCAATAAAATA
>S2:Contig_1
AATGATATAATACTTAATA
>S2:Contig_2
AAAATATTAAAGAAATAATGAAAGAAA
>S3:Contig_1
ATAAAGACGATAAAATAATAAAAGCCAAATCCGACAAAGAAAGAACGG
>S3:Contig_2
AATATTTTAGAGAAAGACATAAACAATAAGAAAAGTATT
>S3:Contig_3
CAAAT
```
### Easy multi binning mode
```bash
SemiBin multi_easy_bin -i combined.fna -b S1.bam S2.bam S3.bam S4.bam S5.bam -o multi_output
```

### Advanced workflows

(1)  Generate `data.csv/data_split.csv` 
```bash
SemiBin generate_data_multi -i combined.fna -b S1.bam S2.bam S3.bam S4.bam S5.bam -o output -s :
```
(2) Generate cannot-link 
```bash
SemiBin generate_cannot_links -i S1.fna -o S1_output
```
```bash
SemiBin generate_cannot_links -i S2.fna -o S2_output
```
```bash
SemiBin generate_cannot_links -i S3.fna -o S3_output
```
```bash
SemiBin generate_cannot_links -i S4.fna -o S4_output
```
```bash
SemiBin generate_cannot_links -i S5.fna -o S5_output
```
(3) Train
```bash
SemiBin train -i S1.fna --data multi_output/samples/S1/train.csv --data-split multi_output/samples/S1/train_split.csv -c S1_output/cannot/cannot.txt -o S1_output --mode single
```
```bash
SemiBin train -i S2.fna --data multi_output/samples/S2/train.csv --data-split multi_output/samples/S2/train_split.csv -c S2_output/cannot/cannot.txt -o S2_output --mode single
```
```bash
SemiBin train -i S3.fna --data multi_output/samples/S3/train.csv --data-split multi_output/samples/S3/train_split.csv -c S3_output/cannot/cannot.txt -o S3_output --mode single
```
```bash
SemiBin train -i S4.fna --data multi_output/samples/S4/train.csv --data-split multi_output/samples/S4/train_split.csv -c S4_output/cannot/cannot.txt -o S4_output --mode single
```
```bash
SemiBin train -i S5.fna --data multi_output/samples/S5/train.csv --data-split multi_output/samples/S5/train_split.csv -c S5_output/cannot/cannot.txt -o S5_output --mode single
```
(4) Bin
```bash
SemiBin bin -i S1.fna --model S1_output/model.h5 --data multi_output/samples/S1/data.csv -o output 
```
```bash
SemiBin bin -i S2.fna --model S2_output/model.h5 --data multi_output/samples/S2/data.csv -o output
```
```bash
SemiBin bin -i S3.fna --model S3_output/model.h5 --data multi_output/samples/S3/data.csv -o output 
```
```bash
SemiBin bin -i S4.fna --model S4_output/model.h5 --data multi_output/samples/S4/data.csv -o output 
```
```bash
SemiBin bin -i S5.fna --model S5_output/model.h5 --data multi_output/samples/S5/data.csv -o output 
```

