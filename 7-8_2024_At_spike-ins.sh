#!/bin/bash
#SBATCH --job-name=At_spike-ins                                               # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7-8_2024_At/At_spike-ins.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/7-8_2024_At/At_spike-ins.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7-8_2024_At
mkdir spike_ins
for file in Raw_Data/*R1*
do
    file2="${file:9:-13}"
    zcat "$file" | sed -n '2~4p' | grep TGCAAATAGGCGGCC | sed -n -e 's/AAAAAAAAAA.*/AAAA/p' | sort | uniq | sed -n -e 's/^.*TGCAAATAGGCGGCC//p' | cut -c1-12 | sort | uniq -c | sort -nr | head -n 100 > spike_ins/raw_"$file2".txt
done