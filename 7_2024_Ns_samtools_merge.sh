#!/bin/bash
#SBATCH --job-name=Ns_samtools_merge                                      # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=2:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7_2024_Ns/Ns_samtools_merge.out        # Location of standard output file
#SBATCH --error=/scratch/jms53460/7_2024_Ns/Ns_samtools_merge.err         # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7_2024_Ns
module load SAMtools/1.16.1-GCC-11.3.0
samtools merge -@ 6 -o 7_2024_Ns_merged.bam -b bamlist
samtools sort -@ 6 7_2024_Ns_merged.bam -o 7_2024_Ns_merged_s.bam
samtools index 7_2024_Ns_merged_s.bam
