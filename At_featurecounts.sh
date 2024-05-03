#!/bin/bash
#SBATCH --job-name=At_featurecounts                                           # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/4_2024_At_Spike_ins/At_featurecounts.out   # Location of standard output file
#SBATCH --error=/scratch/jms53460/4_2024_At_Spike_ins/At_featurecounts.err    # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/4_2024_At_Spike_ins
ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate ./subread-env/

featureCounts -T 6 -s 1 -a TAIR10.1_Col_chr.gff -t 'gene' -g 'ID' -o stringtie_out/read_counts.tab --readExtension5 500 -R BAM hisat2_out/merged_s.bam

featureCounts -T 6 -s 1 -a stringtie_out/sorted_SNPsplit.gtf -o stringtie_out/read_counts.tab --readExtension5 500 -R BAM hisat2_out/merged_s.bam

conda deactivate