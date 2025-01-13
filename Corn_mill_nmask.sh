#!/bin/bash
#SBATCH --job-name=Zm_nmask                                               # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Corn_mill/Zm_nmask.out                # Location of standard output file
#SBATCH --error=/scratch/jms53460/Corn_mill/Zm_nmask.err                 # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Corn_mill
ml BEDTools/2.30.0-GCC-12.2.0
bedtools maskfasta -fi Zm_W22.fa -fo Zm_W22_w_B73_N-masked.fa -bed syri.vcf -fullHeader
