#!/bin/bash
#SBATCH --job-name=Corn_mill_Dedup_txt
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=8                                                  # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Corn_mill/Corn_mill_dedup_txt.out       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Corn_mill/Corn_mill_dedup_txt.err        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Corn_mill

for file in "SNPsplit/"*SNPsplit.bam
do
    file2="${file:9:-13}
umi_tools group -I SNPsplit/"$file2"".SNPsplit.bam" --paired --chimeric-pairs discard --unpaired-reads discard --output-bam -S dedup/"$file2"".Dedup.bam"
        samtools view -f 64 -F 4 dedup/"$file2"".Dedup.bam" | awk -F" " '{split($1, a, "_"); split($21, b, ":"); print $3, $4, $8, $9, substr(a[2],1,6), substr(a[3],1,8), $5, length($10),substr(a[4],1,10), substr(a[4],11,11), substr(b[3],1,2)}' > txtfiles/"$file2"".allele_flagged.txt"
done
