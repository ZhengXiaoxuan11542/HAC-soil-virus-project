#PBS -N host
#PBS -l nodes=1:ppn=64
#PBS -l walltime=999:00:00
#PBS -o host.log
#PBS -e host.err

###########trna##############
mkdir ./07-hostpre/23redo-mag/trna/

tRNAscan-SE -A -o  ./07-hostpre/23redo-mag/trna/trnascan-A \
-m ./07-hostpre/23redo-mag/trna/trnascan-A.stats \
-a ./07-hostpre/23redo-mag/trna/trnascan-A.fasta \
-f ./07-hostpre/23redo-mag/trna/trnascan-A.ss \
--thread 10 \
./03-vdetection/4.vibrant_redo/votu/votu.fa

tRNAscan-SE -B -o  ./07-hostpre/23redo-mag/trna/trnascan-B \
-m ./07-hostpre/23redo-mag/trna/trnascan-B.stats \
-a ./07-hostpre/23redo-mag/trna/trnascan-B.fasta \
-f ./07-hostpre/23redo-mag/trna/trnascan-B.ss \
--thread 10 \
./03-vdetection/4.vibrant_redo/votu/votu.fa

cat ./07-hostpre/23redo-mag/trna/trnascan-A.fasta ./07-hostpre/23redo-mag/trna/trnascan-B.fasta > ./07-hostpre/23redo-mag/trna/trnascan.fasta
makeblastdb -in ./07-hostpre/23redo-mag/trna/trnascan.fasta \
-input_type fasta -dbtype nucl \
-title ./07-hostpre/23redo-mag/trna/vir_trna_db \
-out  ./07-hostpre/23redo-mag/trna/vir_trna_db


 for i in $(ls ./10.bin/mag_50_10_drep/dereplicated_genomes/*.fa); 
  do base_name=$(basename ${i})
  blastn -db ./07-hostpre/23redo-mag/trna/vir_trna_db -query ${i} -out ./07-hostpre/23redo-mag/trna/v_mag_match.out/${base_name} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp" -num_threads 10 -evalue 1e-5 -perc_identity 100;
 done
 #输出到一个文件中，以文件名（bin名）为第一列
 awk '{print FILENAME "\t" $0}' *_bin*.fa >> output_file.txt

###########crispr###############
####cctyper以识别crispr
conda activate cctyper
 for i in $(ls ./10.bin/mag_50_10_drep/dereplicated_genomes/*.fa); 
 do  base_name=$(basename ${i})
 cctyper ${i} ./07-hostpre/23redo-mag/crispr/${base_name} -t 28 --prodigal meta ; done
conda deactivate

#有单个的噬菌体序列文件，优先下面这个代码
#awk '/^>/{f=++d".fasta"} {print > f}' votu.fa  > virusdir/*.fa


####spacephaer以匹配####
spacepharer createsetdb ./03-vdetection/4.vibrant_redo/votu/virusdir/*.fasta  ./07-hostpre/23redo-mag/crispr/pharer/targetSetDB ./07-hostpre/23redo-mag/crispr/pharer/tmpFolder
spacepharer createsetdb ./03-vdetection/4.vibrant_redo/votu/virusdir/*.fasta  ./07-hostpre/23redo-mag/crispr/pharer/targetSetDB_rev ./07-hostpre/23redo-mag/crispr/pharer/tmpFolder --reverse-fragments 1

 for i in $(ls ./07-hostpre/23redo-mag/crispr/*/spacers/*.fa);
 do  base_name=$(basename ${i})
 spacepharer easy-predict ${i}  ./07-hostpre/23redo-mag/crispr/pharer/targetSetDB ./07-hostpre/23redo-mag/crispr/result/results_${base_name}.tsv ./07-hostpre/23redo-mag/crispr/pharer/tmpFolder; done


################blastn##############
  for i in $(ls ./10.bin/mag_50_10_drep/dereplicated_genomes/*.fa); 
  do base_name=$(basename ${i})
  blastn -db ./03-vdetection/4.vibrant_redo/votu/votu_db -query ${i} -out ./07-hostpre/23redo-mag/blastn/results_${base_name}.tsv -evalue 1e-3 -perc_identity 70 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp"; done
# records with alignment length < 2500bp and bitscore < 50 were manually deleted.
