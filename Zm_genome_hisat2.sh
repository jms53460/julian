#!/bin/bash
#SBATCH --job-name=Zm_genome_hisat2                     # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Zm_genome_hisat2.out # Location of standard output file
#SBATCH --error=/scratch/jms53460/Zm_genome_hisat2.err  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022
ml HISAT2/3n-20201216-gompi-2022a
hisat2-build B73_v5_genome.fa.gz B73_v5_hisat2_index

module load SAMtools/1.16.1-GCC-11.3.0
hisat2 -p 6 --dta -x B73_v5_hisat2_index -U A188_genome.fa.gz | samtools view -bS -> hisat2_out/A188_B73_unsorted.bam
samtools sort -@ 6 hisat2_out/A188_B73_unsorted.bam -o A188_B73_s.bam