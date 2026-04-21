#!/bin/bash
#SBATCH --job-name=MaizeMu_4_2026_features_UMIs                           # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/MaizeMu_4_2026/featUMItools.out        # Location of standard output file
#SBATCH --error=/scratch/jms53460/MaizeMu_4_2026/featUMItools.err         # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/MaizeMu_4_2026/

mkdir featurecounts
mkdir bams
mkdir UMIcounts
mkdir UMIcounts_g1
mkdir UMIcounts_g2
ml Mamba/23.11.0-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -t 'gene' -g 'ID' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM hisat2_out/*_s.bam 

conda deactivate

ml SAMtools/1.21-GCC-13.3.0
module load UMI-tools/1.1.4-foss-2023a
for file in "featurecounts/"*_s.bam*
do
    file2="${file:14:-16}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done