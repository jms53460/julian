#!/bin/bash
#SBATCH --job-name=Ns_6_24_download                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=1                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=6:00:00                                                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/June2024Seq/Ns_6_24_download.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/June2024Seq/Ns_6_24_download.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cp /work/bnlab/June2024Sequencing/NTD* /scratch/jms53460/June2024Seq