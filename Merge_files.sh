#!/bin/bash
#SBATCH --job-name=combiningfiles
#SBATCH --partition=batch                   # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=8                   # Number of CPU cores per task
#SBATCH --mem=100gb                         # Job memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec
#SBATCH --output=combofiles_output
#SBATCH --error=combofiles_error
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jms53460@uga.edu        # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd /scratch/jms53460/Ran
mkdir Raw_Data_merged
for file in "Raw_Data/"*L004_R1_001.fastq.gz
do
    file2="${file:9:-25}"

cat "Raw_Data/""$file2""*L004_R1_001.fastq.gz" "Raw_Data/""$file2""*L008_R1_001.fastq.gz" > "Raw_Data_merged/""$file2""R1.fastq.gz"
cat "Raw_Data/""$file2""*L004_R2_001.fastq.gz" "Raw_Data/""$file2""*L008_R2_001.fastq.gz" > "Raw_Data_merged/""$file2""R2.fastq.gz"

done