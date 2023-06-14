

```bash
1. library construction
2. BLASPT
1) 
$ diamond makedb --in qs_refer_seqs.fasta -d qsdb
$ diamond blastp -d qsdb.dmnd -q /publicgp/suying/metagenome_oralbiofilm_output/6.0-specie-annotation/protein.fa -e 1e-5 -f 6 -o 2out_diamond.m6 -k 1 -p 30 

2) 
$ makeblastdb -in qs_refer_seqs.fasta -dbtype prot -parse_seqids -out qs -title qs -logfile qs
$ blastp -query /publicgp/suying/metagenome_oralbiofilm_output/6.0-specie-annotation/protein.fa -out oral_qs_enzyme_blast -db refer_seqs/qs -evalue 1e-5 -outfmt "6 std qlen slen" -max_target_seqs 1 -num_threads 38

$ sed '1i\gene\tseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen' oral_qs_enzyme_blast > titled_oral_qs_enzyme_blast

$ vim -b mSystems-QSseq-ID-classification.txt 
$ %s/\r//g 

$ awk -F '\t' 'NR==FNR{a[$6]=$0;next}NR>FNR{if($2 in a) print a[$2]"\t"$0}' refer_seqs/mSystems-QSseq-ID-classification.txt titled_oral_qs_enzyme_blast > qsblast_info.txt 

$ less oral_qs_enzyme_blast |awk {'print NF'}|sort -u 

3）
$ ssh cn1
$ source activate cdnit
$ makeblastdb -in recA_db.fasta -dbtype prot -parse_seqids -out oral_recA_blast -title recA -logfile recA
$ blastp -query /publicgp/suying/metagenome_oralbiofilm_output/6.0-specie-annotation/protein.fa -out oral_recA_blast -db recA -evalue 1e-5 -outfmt "6 std qlen slen" -max_target_seqs 1 -num_threads 48
$ sed '1i\gene\tseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen' oral_qs_enzyme_blast > titled_oral_qs_enzyme_blast
```

```bash
$ awk -F '\t' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a) print $0}' 1_AHK_rep-14 final-rela0abun-anno.txt > 1_AHK_rep-14_info.txt
```

```bash
$ awk 'BEGIN{OFS=FS="\t"}{if($0~/>/) name=$0 ;else seq[name]=seq[name]$0;}END{for(i in seq) {if(length(seq[i])>0) print i"\n"seq[i]}}' protein.fa > protein_0.fa
$ less protein_0.fa |grep '>'|cut -d '>' -f 2 > protein_lth_morethan_0.txt 
```

