

# conda activate atac-seq
# nohup bash bam_to_bigwig.sh &> out/bam_to_bw.out &

# the reference file path, change to where you put the reference
WDIR="/home/psurana/projects/Matei_lab/chip_cut_run/"
REF="/home/psurana/projects/Matei_lab/bowtie2/hg38"
RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res"
DATA_DIR="/home/psurana/projects/Matei_lab/chip_cut_run/data/raw/"



cd ${RES}/bams/dups_removed/

## bam to bigwig
for i in *.bam; 
do
	bamCoverage -b ${i} -o  ${RES}/bigwig/${i}.bw
done