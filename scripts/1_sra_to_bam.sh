#!/bin/bash
set -euo pipefail

# Initialize directories and reference paths
WDIR="/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/"
REF="/home/psurana/projects/Matei_lab/bowtie2/hg38"
RES="${WDIR}/res"
DATA_DIR="${WDIR}/raw/fastq"

# Function to run FastQC and MultiQC
run_qc() {
    local dir=$1
    local out_dir=$2
    local qc_report=$3
    fastqc "$dir"/*.fastq.gz -o "$out_dir"
    cd "$out_dir"
    export LC_ALL=en_US.utf-8
    export LANG=en_US.utf-8
    multiqc . -o "$qc_report"
}

# Create directories
mkdir -p ${RES}/{bams,fastqc,motifs,peaks}

# Activate bioinformatics environment
conda activate bioinfo

##-----------------------------------------------------------------------------

# Pre-trim quality control
run_qc "${DATA_DIR}" "${RES}/fastqc" "pre_multiqc.html"

##-----------------------------------------------------------------------------
# Adapter trimming
cd "${DATA_DIR}"
for i in *.fastq.gz; do
    cutadapt --minimum-length 20 -a CTGTCTCTTATACACATCT -o "${RES}/trim_fq/$i" $i
done

# Post-trim quality control
run_qc "${RES}/trim_fq" "${RES}/fastqc/post_trim/" "post_multiqc.html"

##-----------------------------------------------------------------------------

# Read alignment using Bowtie2
cd "${RES}/trim_fq"
for sample in *.fastq.gz; do
    base=$(basename "$sample" "_R1_001.fastq.gz")
    bowtie2 -x ${REF}/hg38 -U ${sample} -S ${base}.sam
done
mv *.sam ${RES}/sams/


##-----------------------------------------------------------------------------

# Function for SAM validation using Picard
validate_sam() {
    cd ${SAMS_DIR}
    for sam_file in OV5-DN-1_S4 OV5-DN-2_S5 OV5-DN-3_S6 OV5-DP-1_S1 OV5-DP-2_S2 OV5-DP-3_S3; do
        picard -Xmx5g ValidateSamFile I=${sam_file}.sam M=SUMMARY O=${sam_file}.txt
    done
}

##-----------------------------------------------------------------------------

# Function for filtering, sorting, and indexing SAM/BAM files
process_sam_bam() {
    cd ${SAMS_DIR}
    
    # Remove unwanted chromosomes and create filtered SAM files
    for sam in *.sam; do
        sed '/chrM/d;/random/d;/chrUn/d;/EBV/d' < ${sam} > ${sam}_filtered.sam
    done
    
    # Sort and convert to BAM, then index
    cd ${SAMS_DIR}/filtered_sam
    for filtered_sam in *.sam_filtered.sam; do
        samtools sort $filtered_sam -o ${BAMS_DIR}/${filtered_sam}.bam
        samtools index ${BAMS_DIR}/${filtered_sam}.bam
    done
    
    # Rename files for cleaner names
    cd ${BAMS_DIR}
    for f in *.sam_filtered.sam.bam; do
        mv -- "$f" "${f%.sam_filtered.sam.bam}.bam"
    done
    for f in *.sam_filtered.sam.bam.bai; do
        mv -- "$f" "${f%.sam_filtered.sam.bam.bai}.bam.bai"
    done
}

##-----------------------------------------------------------------------------

# Function for marking duplicates and calculating statistics
mark_and_stats() {
    cd ${BAMS_DIR}
    for bam in *.bam; do
        picard -Xmx5g MarkDuplicates I=${bam} O=${bam}_duplicates.bam M=${bam}_dup_metrics.txt
        samtools flagstat ${bam}_duplicates.bam > ${BAMS_DIR}/bam_stat/${bam}_duplicates.stat
    done
}

##-----------------------------------------------------------------------------

# Validate SAM files
validate_sam

# Process SAM to filtered, sorted, and indexed BAM
process_sam_bam

# Mark duplicates and calculate stats
mark_and_stats

# End of script

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
