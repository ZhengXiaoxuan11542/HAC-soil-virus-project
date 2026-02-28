#PBS -N vdetect
#PBS -l nodes=1:ppn=64 
#PBS -l walltime=999:00:00
#PBS -o vdetect.log
#PBS -e vdetect.err

#1.vibrant 
source /public/home/xiaoxuan/miniconda3/etc/profile.d/conda.sh
conda activate py37
for ((i=1;i<=9;i++)); do 
VIBRANT_run.py -t 40 -l 5000 \
-i ./02-assembly/202310redo/mixassem.rename.5kbsubseq.fa.split/mixassem.rename.5kbsubseq.part_00${i}.fasta  \
-folder ./03-vdetection/4.vibrant_redo/ \
-d miniconda3/pkgs/vibrant-1.2.1-hdfd78af_2/share/vibrant-1.2.1/db/databases/ ; done

conda deactivate 


#2.vs2
conda activate vs2
mkdir ./03-vdetection/4.vibrant_redo/vs2check
virsorter run --keep-original-seq \
-i ./03-vdetection/4.vibrant_redo/vcontig.fna \
-w ./03-vdetection/4.vibrant_redo/vs2check \
--include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae \
--min-length 5000 --min-score 0.5 -j 40 all
conda deactivate



#3.vibrant re-check
conda activate py37
VIBRANT_run.py -t 40 -virome \
-i ./03-vdetection/4.vibrant_redo/vb-recheck/final-viral-combined.rename.fa  \
-folder ./03-vdetection/4.vibrant_redo/vb-recheck \
-d miniconda3/pkgs/vibrant-1.2.1-hdfd78af_2/share/vibrant-1.2.1/db/databases/
conda deactivate



#4.dereplication
cd-hit -i ./03-vdetection/4.vibrant_redo/vb-recheck/VIBRANT_final-viral-combined.rename/final-viral-combined.rename.phages_combined.fna \
-c 0.95 -aS 0.85 -T 20 -M 35000 \
-o ./03-vdetection/4.vibrant_redo/votu/votu.fa



#5.checkv summary
checkv end_to_end ./03-vdetection/4.vibrant_redo/votu/votu.fa ./03-vdetection/4.vibrant_redo/votu/checkv_out -t 20



#6.viral classification
conda activate genomad
genomad annotate ./03-vdetection/4.vibrant_redo/votu/votu.fa ./04-vannotation/23redo/genomad_out database/genomad_db/genomad_db
conda deactivate



#7. AMG identification	
#7.1VIBRANT
conda activate py37
VIBRANT_run.py -virome -t 50 -i ./03-vdetection/4.vibrant_redo/votu/votu.fa \
-folder ./05.23redo_vamg/vibrant \
-d miniconda3/pkgs/vibrant-1.2.1-hdfd78af_2/share/vibrant-1.2.1/db/databases/
conda deactivate

#7.2DRAM-v
conda activate vs2
virsorter run --prep-for-dramv -w ./05.23redo_vamg/vs2-dmv \
-i ./03-vdetection/4.vibrant_redo/votu/votu.fa \
--min-length 5000 --min-score 0.5 -j 50 all
conda deactivate

conda activate DRAM
DRAM-v.py annotate -i  ./05.23redo_vamg/vs2-dmv/for-dramv/final-viral-combined-for-dramv.fa \
-v ./05.23redo_vamg/vs2-dmv/for-dramv/viral-affi-contigs-for-dramv.tab \
-o ./05.23redo_vamg/dmv \
--skip_trnascan --threads 50 --min_contig_size 5000
DRAM-v.py distill -i ./05.23redo_vamg/dmv/annotations.tsv -o ./05.23redo_vamg/dmv/dramv-distill
conda deactivate

#7.1 and #7.2 output information were sorted together

#8. Coverage table generation
#Build BAM file 
READ_DIR=mydata/hnmeta/01-qc/3.rmhost
reads=$(for r1 in $READ_DIR/*.nohost.R1.fq; do
    echo -n "$r1 ${r1/.R1.fq/.R2.fq} "
done)

coverm make --threads 40 \
--reference ./03-vdetection/4.vibrant_redo/votu/votu.fa \
--coupled $reads \
--output-directory ./05-abundance/redo-abundance/

#Filtered out reads which percent identity with less than 95%
coverm filter -b ./05-abundance/redo-abundance/votu.fa.CK-1-J.nohost.R1.fq.bam -o ./05-abundance/redo-abundance/CK-1-J_filtered.bam  --min-read-percent-identity 0.95 --threads 40

#Abundance table construction
coverm contig --methods trimmed_mean \
-b ./05-abundance/redo-abundance/*_filtered.bam  \
--min-covered-fraction 0.75 >  ./05-abundance/redo-abundance/votu.tmean.tsv



