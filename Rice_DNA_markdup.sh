#!/bin/bash
#SBATCH --job-name=Rice_DNA_markdup                                # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=100gb                                                # Total memory for job
#SBATCH --time=12:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_DNA/Rice_markdup.out       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_DNA/Rice_markdup.err        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_DNA
module load SAMtools/1.16.1-GCC-11.3.0
mkdir markdup

#This first collate command can be omitted if the file is already name ordered or collated:
samtools collate -o markdup/namecollate.bam hisat2_out/NRE1-2_s.bam
#Add ms and MC tags for markdup to use later:
samtools fixmate -m markdup/namecollate.bam markdup/fixmate.bam
#Markdup needs position order:
samtools sort -o markdup/positionsort.bam markdup/fixmate.bam
#Finally mark duplicates:
samtools markdup --duplicate-count markdup/positionsort.bam markdup/NRE1-2_markdup.bam