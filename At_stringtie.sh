#!/bin/bash
#SBATCH --job-name=At_stringtie                                           # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/4_2024_At_Spike_ins/At_stringtie.out   # Location of standard output file
#SBATCH --error=/scratch/jms53460/4_2024_At_Spike_ins/At_stringtie.err    # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/4_2024_At_Spike_ins
ml SAMtools/1.16.1-GCC-11.3.0
samtools sort -@ 6 At_SNPsplit/merged_s.allele_flagged.bam -o At_SNPsplit/sorted_SNPsplit.bam
mkdir stringtie_out

ml StringTie/2.2.1-GCC-11.3.0
stringtie At_SNPsplit/sorted_SNPsplit.bam -p 6 -G TAIR10.1_Col_5.gff --rf -o stringtie_out/sorted_SNPsplit.gtf
