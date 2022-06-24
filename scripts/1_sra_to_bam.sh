#! /bin/bash
set -euo pipefail

# # scripts already stored in scripts dir

#- nohup bash 1_sra_to_bam.sh &> out/6_bam2dups.out &

##- total 23 files

##* make dirs
# mkdir -p results/{bams,fastqc,motifs,peaks}
# #mkdir -p hg38
# mkdir -p input

# conda activate atac-seq

##--------------------------------------------------------------------------

# the reference file path, change to where you put the reference
WDIR="/home/psurana/projects/Matei_lab/chip_cut_run/"
REF="/home/psurana/projects/Matei_lab/bowtie2/hg38"
RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res"
DATA_DIR="/home/psurana/projects/Matei_lab/chip_cut_run/data/raw/"

## extract the fastq files - not needed as fastq files are provided
# fasterq-dump $SRA


# ##---------------------------------------------------------------------------

# #* step1, quality control of the fastq files
# echo "Fastqc quality control usually  takes long -- Be patient..."
# cd ${DATA_DIR}
# fastqc *.fastq.gz -o $RES/fastqc

# # ##---------------------------------------------------------------------------

# ##* step 1a, combine all fastqc reports using multiqc
# export LC_ALL=en_US.utf-8
# export LANG=en_US.utf-8
# cd ${RES}/fastqc/1_pre_trim/ 
# multiqc . -o pre_multiqc.html

##---------------------------------------------------------------------------

# ##* step 1b, Adapter trimming - Nextera transposase here
# # took over an hour here i think
# cd ${DATA_DIR}
# for i in *.fastq.gz; 
# do 
# cutadapt --minimum-length 20 -a CTGTCTCTTATACACATCT -o ${RES}/trim_fq/${i} $i
# done
# echo "Trimming of adapters done by cutadapt"

# ##* step 1c, run fastqc again and see if adapters have gone!
# cd ${RES}/trim_fq
# fastqc *.fastq.gz -o $RES/fastqc/2_post_trim/

# export LC_ALL=en_US.utf-8
# export LANG=en_US.utf-8

# cd ${RES}/fastqc/2_post_trim/
# multiqc . -o post_multiqc.html

# ##---------------------------------------------------------------------------

# ##* step2, align the reads with bowtie2

# for sample in `ls /home/psurana/projects/Matei_lab/chip_cut_run/data/res/trim_fq/*.fastq.gz`
# do
# dir="/home/psurana/projects/Matei_lab/chip_cut_run/data/res/trim_fq"
# res="/home/psurana/projects/Matei_lab/chip_cut_run/data/res/sams"
# base=$(basename $sample "_001.fastq.gz")
# bowtie2 -x ${REF}/hg38 -U ${dir}/${base}_001.fastq.gz -S ${res}/${base}.sam --no-unal
# done
# mv *.sam ${RES}/sams/

##-----------------------------------------------------------------------------

# conda activate rna-seq

# java -jar ./anaconda3/pkgs/picard-2.18.29-0/share/picard-2.18.29-0/picard.jar

#- validate sam file
# cd ${RES}/sams/

# f=(*sam)
# echo "${f[0]}"

# for i in {0..21}
# do
# picard -Xmx5g ValidateSamFile I=${f[i]} M=SUMMARY O=${f[i]}.txt
# done

# picard -Xmx5g ValidateSamFile I=${f[0]} M=VERBOSE O=${f[0]}.txt

#------------------------------------------------------------------------------

# ##* step 3a, convert sam to bam, and index the bam

#------------------------------------------------------------------------------
## ref
## sed '/chrM/d;/random/d;/chrUn/d;/EBV/d' < OV5-DN-1_S4.sam  > filter_OV5-DN-1_S4.sam 
## samtools sort filter_OV5-DN-1_S4.sam -o filter_OV5-DN-1_S4.bam
## samtools index filter_OV5-DN-1_S4.bam
## #- check
## samtools idxstats filter_OV5-DN-1_S4.bam| cut -f 1 
#------------------------------------------------------------------------------

# ##* step 3a, convert sam to bam, and index the bam

#--- remove other chromosomes
# cd ${RES}/sams/
# for i in *.sam; 
# do
# sed '/chrM/d;/random/d;/chrUn/d;/EBV/d' < ${i} > ${i}_filtered.sam
# done


# mkdir filtered_sam
# mv *.sam_filtered.sam filtered_sam/

# #--- filter sam file 
# cd ${RES}/sams/filtered_sam
# for i in *.sam_filtered.sam; 
# do
# samtools sort $i -o ${RES}/bams/$i.bam
# done


# cd ${RES}/bams/
# for i in *.bam;
# do
# samtools index $i
# done

# #- rename file extensions
# for f in *.sam_filtered.sam.bam; 
# do 
#     mv -- "$f" "${f%.sam_filtered.sam.bam}.bam"
# done
# for f in *.sam_filtered.sam.bam.bai; 
# do 
#     mv -- "$f" "${f%.sam_filtered.sam.bam.bai}.bam.bai"
# done


# ##---------------------------------------------------------------------------

#- 3b, remove duplciates in chip seq datasets, calculate and print statistics of alignment
#- calculate the dups

# # conda activate rna-seq
# cd ${RES}/bams/
# for i in *.bam ; do 
#     picard -Xmx5g MarkDuplicates I=${i} O=${RES}/bams/dups_removed/${i}_duplicates.bam M=${i}_dup_metrics.txt REMOVE_SEQUENCING_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
# done


# ##---------------------------------------------------------------------------


#- samtools - alignment stats to see if the bam file should be cleaned further

# for f in *.bam_duplicates.bam; 
# do 
#     mv -- "$f" "${f%.bam_duplicates.bam}_duplicates.bam"
# done

# cd ${RES}/bams/dups_removed
# for i in *.bam ; do 
#     samtools flagstat "$i" > ${RES}/bams/dup_bam_stat/"$i".stat
# done

# cd ${RES}/bams/dups_removed
# for i in *.bam;
# do
# samtools index $i
# done

# #-overall looks good
# ##---------------------------------------------------------------------------

# ## ----- end of script ------
