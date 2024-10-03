#!/bin/bash
#SBATCH --job-name=Ns_features_UMIs                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7_2024_Ns/Ns_features_UMIs3.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/7_2024_Ns/Ns_features_UMIs3.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7_2024_Ns
mkdir featurecounts3
ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate ./subread-env/

featureCounts -T 6 -s 2 -a nsyl_10.gtf -o featurecounts3/read_counts.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit.bam

conda deactivate