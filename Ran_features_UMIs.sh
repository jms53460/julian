#!/bin/bash
#SBATCH --job-name=Ran_features_UMIs                                     # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Ran/Ran_features_UMIs.out              # Location of standard output file
#SBATCH --error=/scratch/jms53460/Ran/Ran_features_UMIs.err               # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Ran

mkdir featurecounts
mkdir bams
mkdir UMIcounts
ml Mamba/23.11.0-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a JEC21_56.gff -t 'protein_coding_gene' -g 'gene_id' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM Mapped/*.bam

conda deactivate

ml SAMtools/1.21-GCC-13.3.0
module load UMI-tools/1.1.4-foss-2023a
for file in "featurecounts/".bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done