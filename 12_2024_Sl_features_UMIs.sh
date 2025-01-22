#!/bin/bash
#SBATCH --job-name=Features_UMIs                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=70gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/12_2024_Sl/Features_UMIs.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/12_2024_Sl/Features_UMIs.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/12_2024_Sl
for file in "featurecounts2/"*g2.bam*
do
    file2="${file:15:-22}"
    if [ ! -f "UMIcounts2_g2/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams2/$file2"
        samtools index "bams2/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams2/$file2" -S "UMIcounts2_g2/${file2}.tsv"
    fi
done