#!/bin/bash
#SBATCH --job-name=Ns_features                          # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_features.out         # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_features.err          # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460
ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate /home/jms53460/./subread-env/

featureCounts -T 6 -s 1 -a "stringtie_out2/stringtie_merged.gtf" -o "stringtie_out2/read_counts.tab" --readExtension5 500 -R BAM "hisat2_out/"*.bam

conda deactivate
