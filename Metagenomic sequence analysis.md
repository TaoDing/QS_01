

```bash
adapter：
R1-adapter:	GATCGGAAGAGCACACGTCTGAACTCCAGTC
R2-adapter:	GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```

## 1. Data quality control

### 1.1 qc Quality 

```bash
1.qc Quality 
[suying@mgt fq]$ ls *gz>ID.txt
[suying@mgt fq]$ for i in `cat ID.txt`;do fastqc ${i} -o fastqc -t 12;done
[suying@mgt fq]$ nohup sh -c 'for i in `cat ID.txt`;do fastqc ${i} -o fastqc -t 12;done' &
```

```bash
awk -F. '{print $1}' ID.txt > ID.new.txt 
$ conda create -n Cutadapt python=3.7 
$ conda activate Cutadapt  
$ conda deactivate  
$ conda install -c bioconda cutadapt

$ for i in `cat ID.new3.txt`;do cutadapt  -a GATCGGAAGAGCACACGTCTGAACTCCAGTC -A GATCGGAAGAG                   CGTCGTGTAGGGAAAGAGTGT -e 0.1 -O 10 -m 100 -n 1 -j 8 --pair-filter=both -o ${i}.noadapter.R1.fastq -p ${i}.noadapter.R2.fastq   ${i}.R1.fq.gz ${i}.R2.fq.gz;done
```

```bash
3.
$ conda create -n kneaddata
$ conda activate kneaddata
$ conda install -c biobakery kneaddata

$ mkdir -p db
$ kneaddata_database --download human_genome bowtie2 db/
$ tar -zxvf Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz

parallel -j 3 --xapply 'echo {1} {2}'    ::: clean.data/*1.clean.fq.gz ::: clean.data/*2.clean.fq.gz
   $ nohup parallel -j 2 --xapply 'kneaddata -i {1} -i {2} -o kneaddata_out2 -v -db /public7/suying/metagenome_oralbiofilm_output/db  --trimmomatic /public/apps/anaconda2/bin/Trimmomatic-0.38/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" -t 24 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output'  

 $ nohup sh -c 'for i in `cat knID-new.txt`;do kneaddata -i ${i}.noadapter.R1.fastq -i ${i}.noadapter.R2.fastq -o kneaddata_out2 -v -db /public7/suying/metagenome_oralbiofilm_output/db  --trimmomatic /public/apps/anaconda2/bin/Trimmomatic-0.38/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" -t 24 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output;done' &
```

## 2. Megahit 

```bash
$ conda create -n megahit python=3.7
$ conda install -c bioconda megahit
```

```bash

$ nohup megahit -1 D10J1.noadapter.R1_kneaddata_paired_1.fastq  -2 D10J1.noadapter.R1_kneaddata_paired_2.fastq  -o test2/ --k-list 21,29,39,59,79,99,119,141 -m 0.9 -t 20 & 

$ ls *.fastq>name.txt
$ awk -F. '{print $1}' name.txt |uniq>nameid.txt
$ nohup sh -c 'for i in `cat nameid.txt`;do megahit -1 ${i}.noadapter.R1_kneaddata_paired_1.fastq  -2 ${i}.noadapter.R1_kneaddata_paired_2.fastq  -o assemble1/ --k-list 21,29,39,59,79,99,119,141 -m 0.9 -t 20;done' &

```

```bash
6. QUAST
$ quast.py result/megahit/final.contigs.fa -o result/megahit/quast -t 12 
$ time metaquast.py result/megahit/final.contigs.fa -o result/megahit/metaquast -t 12
```

## 3. Uniqe gene

### 3.1 metaProdigal 

```bash

 $ mkdir -p temp/prodigal 
 $ nohup prodigal -i assemble1/final.contigs.fa -d temp/prodigal/D1_gene.fa -o temp/prodigal/D1_gene.gff -p meta -f gff >temp/prodigal/D1_gene.log 2>&1 & 
 $ tail temp/prodigal/gene.log 
 $ seqkit stats temp/prodigal/gene.fa 
 $ grep -c 'partial=00' temp/prodigal/gene.fa 
 
 $ grep 'partial=00' temp/prodigal/gene.fa | cut -f1 -d ' '| sed 's/>//' > temp/prodigal/full_length.id
 $ seqkit grep -f temp/prodigal/full_length.id temp/prodigal/gene.fa > temp/prodigal/full_length.fa
 $ seqkit stat temp/prodigal/full_length.fa 
```

### 3.2 cd-hit

```bash
$ conda create -n cdhit
$ source activate cdhit

$ mkdir -p result/NR

$ cp D1_gene.fa ./ 

$ sed -i 's/>/>D1_/' ./D1_gene.fa
$ sed 's/ .*//' D1_gene.fa > D1_gene_newname.fa

$ cat *_newname.fa > all.fa

$ nohup cd-hit-est -i allgene.fa -o result/NR/nucleotide.fa -aS 0.9 -c 0.95 -G 0 -g 1 -T 0 -M 0 -d 0 &

$ grep -c '>' result/NR/nucleotide.fa
$ seqkit translate --trim result/NR/nucleotide.fa > result/NR/protein.fa  

$ awk 'BEGIN{OFS=FS="\t"}{if($0~/>/) name=$0 ;else seq[name]=seq[name]$0;}END{for(i in seq) {if(length(seq[i])>100) print i"\n"seq[i]}}' protein.fa > protein_seq_100.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < nucleotide_seq.fasta > nucleotide_seq_new.fasta
grep -A1 -Fw -f protein_seq_102.fasta nucleotide_seq_new.fasta | grep -v -- "--" > nucleotide_seq_new_102.fasta  
```

### 3.3 salmon

```bash
1. 
$ conda activate cdhit
$ bwa index /public7/suying/metagenome_oralbiofilm_output/4.0-genecdhit/result/NR/nucleotide.fa 
$ ls *fastq>name.txt
$ for i in `cat name.txt`;do fastq-sort $i > sort_$i ;done 
or
$ fastq-sort D10J1.noadapter.R1_kneaddata_paired_1.fastq > test_D10J1.noadapter.R1_kneaddata_paired_1.fastq 
$ fastq-sort D10J1.noadapter.R1_kneaddata_paired_2.fastq > test_D10J1.noadapter.R1_kneaddata_paired_2.fastq
$ grep '@' D10J1.noadapter.R1_kneaddata_paired_1.fastq|tail -n 10 
$ grep '@' D10J1.noadapter.R1_kneaddata_paired_2.fastq|tail -n 10 

$ ls *fastq|awk -F. '{print $1}'|uniq>ID.txt
$ nohup sh -c 'for i in `cat ID.txt`; do /public/home/suying/.conda/envs/cdhit/bin/bwa mem -t 24 /public7/suying/metagenome_oralbiofilm_output/4.0-genecdhit/result/NR/nucleotide.fa sort_$i.noadapter.R1_kneaddata_paired_1.fastq sort_$i.noadapter.R1_kneaddata_paired_2.fastq> sort_$i.aln.sam;done' & 

```

```bash

3.2 
$ for i in *.aln.sam; do grep -v ^@ $i | cut -f 3 | sort |uniq -c| awk 'BEGIN{ FS=" ";OFS="," }{print $2,$1}'| awk 'BEGIN{ FS=",";OFS=","}{ if ($2 > 2) print $1,$2;else print $1,"0"}' > count-out/$i.count.csv;done 

$ grep -v ^@ test-D10J1.aln.sam | cut -f 3 | sort |uniq -c| awk 'BEGIN{ FS=" ";OFS="," }{print $2,$1}' | awk 'BEGIN{ FS=",";OFS=","}{ if ($2 > 2) print $1,$2;else print $1,"0"}' > test-D10J1.count.csv 

3.3 ref : 
$ conda install -c bioconda bioawk
$ cd 4.0-genecdhit/
$ bioawk -c fastx '{print $name, length($seq)}' nucleotide.fa | tr "\t" "," > gene.length.csv  


$ python jiyinfengdu.py $i.count.csv $i > $i.tmp  

$ ls *csv>ID.txt
$ awk -F. '{print $1}' ID.txt |awk -F_ '{print $2}' >name.txt
$ mv name.txt /publicgp/suying/metagenome_oralbiofilm_output/5.0-genecaculate/
$ for i in `cat name.txt`; do python jiyinfengdu.py sort_$i.aln.sam.count.csv $i > ${i}.tmp;done
$ paste *.tmp > end.tsv

3.
$ for i in *.tmp;do sed 's/,/\t/g' ${i}  |cut -f 2 > ${i}_abundance.txt;done 
$ sed 's/,/\t/g' D1J1.tmp  |cut -f 1 >gene-name.txt
$ paste D1J?.*txt D2J?.*txt D3J?.*txt D4J?.*txt D5J?.*txt D6J?.*txt D7J?.*txt D8J?.*txt D9J?.*txt D10J?.*txt D11J?.*txt > end.txt
$ paste gene-name.txt end.txt >final_relative_abundance.txt

$ cd /publicgp/suying/metagenome_oralbiofilm_output/6.0-specie-annotation/
$ cp micro-gene-annotation.txt /publicgp/suying/metagenome_oralbiofilm_output/5.0-genecaculate/count-out-relative_abundance/relative-abundance/
$ sed '1i\gene\tAccession\tTaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' micro-gene-annotation.txt>test.txt  
$ mv test.txt micro-gene-annotation.txt
$ awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a)print a[$1],$0}' final_relative_abundance.txt micro-gene-annotation.txt > final-gene-abundance-tax.txt
$ sed 's/,/\t/g' final-gene-abundance-tax.txt > test.txt
$ mv test.txt final-gene-abundance-tax.txt
```

```bash

$ samtools view -bS --threads 20 contig.sam > contig.bam

$ samtools sort contig.bam -o contig.reads.sorted.bam --threads 20

$ samtools index contig.reads.sorted.bam

$ mkdir bins
$ cp contigs.fasta bins
$ checkm coverage -x fasta -m 20 -t 20 bins contigs_coverage.out contig.reads.sorted.bam
$ checkm coverage -x fasta -m 20 -t 64 checkm_bins contigs_coverage.out checkm_bins/1A_contig_sorted.bam
$ checkm coverage -x fa -m 20 -t 64 bins contigs_coverage.out bins/test-D10J1.aln.sam.bam.sorted.bam

```

## 5. taxonomy

```bash
# Generate report in default taxid output

$ DBNAME=/publicgp/suying/references/kraken2_database/db2/
$ mkdir -p DBNAME
$ cd $DBNAME
 
$ kraken2-build  --threads 56 --download-taxonomy --db $DBNAME 

kraken2-build  --threads 56  --download-library bacteria --db $DBNAME
kraken2-build  --threads 56  --download-library viral --db $DBNAME
kraken2-build  --threads 56  --download-library archaea --db $DBNAME
kraken2-build  --threads 56  --download-library plasmid --db $DBNAME
kraken2-build  --threads 56  --download-library human --db $DBNAME
kraken2-build  --threads 56  --download-library fungi --db $DBNAME

$ for i in archaea bacteria plasmid viral human fungi; do kraken2-build --download-library $i --threads 56 --db $DBNAME;done

   
```

```bash

$ conda activate cdhit
$ cp /public7/suying/metagenome_oralbiofilm_output/4.0-genecdhit/result/NR/nucleotide.fa ./
$ mkdir temp 
$ kraken2 --db /public5/pulic4_backups/lym/kracken_database/ss_lym20201231/ \
      nucleotide.fa \
      --threads 12 \
      --report temp/NRgene.report \
      --output temp/NRgene.output
      
# Genes & taxid list
$ grep '^C' temp/NRgene.output|cut -f 2,3|sed '1 i Name\ttaxid' > temp/NRgene.taxid

# Add taxonomy
$ awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $1,a[$2]}' \
      taxonomy.txt \
      temp/NRgene.taxid \
      > temp/nucleotide.tax
      
$ memusg -t /conda2/envs/humann3/bin/python3 /db/EasyMicrobiome/script/summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m result/NR/nucleotide.tax \
      -c '2,3,4,5,6,7,8,9' -s ',+,+,+,+,+,+,+,' -n raw \
      -o result/NR/tax
$ wc -l result/NR/tax*|sort -n



$ /usr/bin/python3 summarizeAbundance.py -i test-D10J1.tmp -m temp/nucleotide.tax -c '2,3,4,5,6,7,8,9' -s ',+,+,+,+,+,+,+,' -n raw -o temp/tax
```

