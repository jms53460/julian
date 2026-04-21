#!/bin/bash
#SBATCH --job-name=At_4_2026_features_UMIs                               # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/At_4_2026/features_UMIs.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/At_4_2026/features_UMIs.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/At_4_2026/

ml SAMtools/1.21-GCC-13.3.0
ml UMI-tools/1.1.4-foss-2023a
for file in "featurecounts/"*SNPsplit.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g1.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g1/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g1/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g2.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g2/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g2/${file2}.tsv"
    fi
done