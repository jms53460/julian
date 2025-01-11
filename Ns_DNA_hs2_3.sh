#!/bin/bash
#SBATCH --job-name=Ns_hisat2                                      # Job name
#SBATCH --partition=batch                                             # Partition (queue) name
#SBATCH --ntasks=1                                                    # Single task job
#SBATCH --cpus-per-task=6                                             # Number of cores per task
#SBATCH --mem=50gb                                                    # Total memory for job
#SBATCH --time=12:00:00                                               # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_h2_3.out       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_h2_3.err        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                  # Where to send mail
#SBATCH --mail-type=END,FAIL                                          # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
mkdir hisat2_out3
cp hisat2_out/TW136xTW137*fastq.gz hisat2_out3

ml HISAT2/3n-20201216-gompi-2022a
hisat2-build TW136_consensus_N-masked_12.fa TW136_consensus_N-masked_12_index
module load SAMtools/1.16.1-GCC-11.3.0
hisat2 -p 6 --dta -x TW136_consensus_N-masked_12_index -1 "hisat2_out3/TW136xTW137_R1_001.fastq.gz" -2 "hisat2_out3/TW136xTW137_R2_001.fastq.gz" | samtools view -bS -> "hisat2_out3/TW136xTW137_unsorted.bam"
samtools sort -@ 6 "hisat2_out3/TW136xTW137_unsorted.bam" -o "hisat2_out3/TW136xTW137_s.bam"
