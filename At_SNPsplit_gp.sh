#!/bin/bash
#SBATCH --job-name=At_SNPsplit_genome_prep                                # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=6:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/4_2024_At_Spike_ins/At_SNPsplit_gp.out # Location of standard output file
#SBATCH --error=/scratch/jms53460/4_2024_At_Spike_ins/At_SNPsplit_gp.err  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/4_2024_At_Spike_ins
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1

SNPsplit_genome_preparation --vcf_file At_vcf.gz --reference_genome Col_genome_dir --strain Ler_0 --skip_filtering --genome_build Col
