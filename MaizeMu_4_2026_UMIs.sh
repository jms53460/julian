#!/bin/bash
#SBATCH --job-name=MaizeMu_4_2026_UMIs                                    # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=200gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/MaizeMu_4_2026/UMIs.out        # Location of standard output file
#SBATCH --error=/scratch/jms53460/MaizeMu_4_2026/UMIs.err         # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/MaizeMu_4_2026/

module load UMI-tools/1.1.4-foss-2023a
for file in "featurecounts/"*_s.bam*
do
    file2="${file:14:-16}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done