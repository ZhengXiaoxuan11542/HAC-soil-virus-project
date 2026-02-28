#PBS -N mag
#PBS -l nodes=1:ppn=60 
#PBS -l walltime=999:00:00
#PBS -o mag.log
#PBS -e mag.err

#1.binning
conda activate metawrap-env
metawrap binning \
-t 60 -m 400 -l 1500 \
--concoct --metabat2  --maxbin2 \
-a ./02-assembly/CK-J-mix/CK-J-mix.1500.relabeld.fa \
-o ./10.bin/CKJ  \
./01-qc/3.rmhost/CKJ.nohost_1.fastq ./01-qc/3.rmhost/CKJ.nohost_2.fastq --run-checkm

#2.bin_refinement
metawrap bin_refinement \
-o ./10.bin/CKJ/binrefin_50_10 \
-A ./10.bin/CKJ/metabat2_bins/ -B ./10.bin/CKJ/maxbin2_bins/ -C ./10.bin/CKJ/concoct_bins/ \
-t 50 -c 50 -x 10
conda deactivate 

#3.rename genome, run in the current folder
prefix="CKJ_"; for file in *; do [ -f "$file" ] && mv "$file" "${prefix}${file}" && echo "已将文件 $file 重命名为 ${prefix}${file}"; done
#sorted together
cp ./10.bin/CKJ/binrefin_50_10/metawrap_50_10_bins/*.fa ./10.bin/bin_50_10/

#4.dereplicate
dRep dereplicate \
./10.bin/mag_50_10_drep/ \
-g ./10.bin/bin_50_10/*.fa \
--S_algorithm ANImf -sa 0.99 -nc 0.1 -comp 50 -con 10 \
--checkM_method lineage_wf --clusterAlg single -p 64



#5.gtdbtk classification
GTDBTK_DATA_PATH=/public/home/xiaoxuan/database/gtdbtk/release207_v2/
gtdbtk classify_wf --genome_dir ./10.bin/mag_50_10_drep/dereplicated_genomes/ \
--out_dir ./10.bin/classify_gtdbtk_output-v207 \
-x .fa --prefix F --cpus 4 --debug \
--skip_ani_screen --pplacer_cpus 1