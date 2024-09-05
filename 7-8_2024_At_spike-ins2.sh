#!/bin/bash
#SBATCH --job-name=At_spike-ins2                                              # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7-8_2024_At/At_spike-ins2.out              # Location of standard output file
#SBATCH --error=/scratch/jms53460/7-8_2024_At/At_spike-ins2.err               # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7-8_2024_At
mkdir spike_ins2
for file in "Demultiplexed/"*s.fastq*
do
    zcat "$file" | sed -n '2~4p' | grep TGCAAATAGGCGGCC | sed -n -e 's/AAAAAAAAAA.*/AAAA/p' | sort | uniq | sed -n -e 's/^.*TGCAAATAGGCGGCC//p' | cut -c1-12 | grep -F -f /scratch/jms53460/7-8_2024_At/96spike_in_barcodes.txt | sort | uniq -c | sort -nr | head -n 22 > spike_ins2/"${file:14:-9}".txt
done

mkdir spike_ins3
for file in spike_ins2/*s.txt
do
    awk '{print $1,$2}' $file OFS="" > spike_ins3/${file:10:-4}.tsv
done
