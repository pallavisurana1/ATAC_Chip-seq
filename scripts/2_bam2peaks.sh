#! /bin/bash
set -euo pipefail

#- conda activate atac-seq

#- nohup bash 2_bam2peaks.sh &> out/peaks_macs2.out &

#- conda install -c bioconda macs2
#- ref - https://biohpc.cornell.edu/doc/epigenomics_2020_exercise2.pdf

WDIR="/home/psurana/projects/Matei_lab/chip_cut_run/"
REF="/home/psurana/projects/Matei_lab/bowtie2/hg38"
RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res"
DATA_DIR="/home/psurana/projects/Matei_lab/chip_cut_run/data/raw/"

##---* step 4 - call peaks with macs2

cd ${RES}/bams/dups_removed
for i in *.bam ; do 
    macs2 callpeak -t ${i} -f BAM -g 2.7e9 -n "${i}" -B -q 0.05 --nomodel --extsize 100 --keep-dup all --call-summits --outdir ${RES}/peaks/macs2
done


# for f in *.bam_peaks.narrowPeak; 
# do 
#     mv -- "$f" "${f%.bam_peaks.narrowPeak}.narrowPeak"
# done


##---* OR ------------
##---* step 4 - call peaks with genrich

#- sort bam files by query name
# cd ${RES}/bams/
# for i in *.bam; 
# do
# samtools sort -n ${i} -o ${RES}/bams/sorted_query_name/${i}
# done


#-genrich
# cd ${RES}/bams/sorted_query_name/
# for i in *.bam ; do 
#     Genrich  -t ${i}  -o $RES/peaks/genrich_p=0.1/narrowPeak/${i}.narrowPeak  -b $RES/peaks/genrich_p=0.1/bed/${i}.bed  \
#     -f $RES/peaks/genrich_p=0.1/${i}.log  -j  -y  -d 150  -p 0.1  -v
# done

# cd $RES/peaks/genrich_p=0.1/narrowPeak/
# for f in *.bam.narrowPeak; 
# do 
#     mv -- "$f" "${f%.bam.narrowPeak}.narrowPeak"
# done

#---*#*----*#*----*#*----*#*----

