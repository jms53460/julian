#!/bin/bash
#SBATCH --job-name=Features_UMIs                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=6:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/11_2024_Sl/Features_UMIs.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/11_2024_Sl/Features_UMIs.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/11_2024_Sl
mkdir UMIcounts3
module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4

for file in "featurecounts2/"*SNPsplit.bam*
do
    file2="${file:15:-22}"
    if [ ! -f "UMIcounts3/${file2}.tsv" ]; then

        umi_tools count -I "bams2/$file2" -S "UMIcounts3/${file2}.tsv"
    
    fi
done