#!/bin/bash
#SBATCH --job-name=At_raw_spike_ins                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=1                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=6:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/4_2024_At_Spike_ins/At_raw_spike_ins.out   # Location of standard output file
#SBATCH --error=/scratch/jms53460/4_2024_At_Spike_ins/At_raw_spike_ins.err    # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/4_2024_At_Spike_ins/Raw_Data
mkdir /scratch/jms53460/4_2024_At_Spike_ins/spike_ins
for file in *R1*
do
    zcat "$file" | sed -n '2~4p' | grep TGCAAATAGGCGGCC | sed -n -e 's/AAAAAAAAAA.+/AAAA/p' | sort | uniq | sed -n -e 's/^.*TGCAAATAGGCGGCC//p' | cut -c1-12 | sort | uniq -c | sort -nr | head -n 100 > /scratch/jms53460/4_2024_At_Spike_ins/spike_ins/raw_"$file".tab
done