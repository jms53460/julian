#!/bin/bash
#SBATCH --job-name=download_7-8_2024_data                                 # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=48:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7-8_2024_At/download.out               # Location of standard output file
#SBATCH --error=/scratch/jms53460/7-8_2024_At/download.err                # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /work/bnlab/July2024Sequencing
sftp somers_10237@dnaseq2.igsp.duke.edu
SMbSq52Nupbq
get -r Somers_10237_240830A9

cd /work/bnlab/Aug2024Sequencing
sftp scroggs_10318@dnaseq2.igsp.duke.edu
t7Chu1Nish49
get -r Scroggs_10318_240830A9

###This did not work. I instead ran these lines in the cluster terminal. (The text between the sftp and get lines are the passwords)
