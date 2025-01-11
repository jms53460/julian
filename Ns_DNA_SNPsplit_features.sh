#!/bin/bash
#SBATCH --job-name=SNPsplit_features                                                # Job name
#SBATCH --partition=batch                                                           # Partition (queue) name
#SBATCH --ntasks=1                                                                  # Single task job
#SBATCH --cpus-per-task=6                                                           # Number of cores per task
#SBATCH --mem=100gb                                                                  # Total memory for job
#SBATCH --time=12:00:00                                                              # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/SNPsplit_features.out                       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/SNPsplit_features.err                        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                                # Where to send mail
#SBATCH --mail-type=END,FAIL                                                        # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
mkdir SNPsplit
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
SNPsplit --conflicting -o SNPsplit --snp_file Ns_SNPs.tab hisat2_out2/TW136xTW137_s.bam
    samtools sort -@ 6 SNPsplit/TW136xTW137_s.allele_flagged.bam -o SNPsplit/TW136xTW137_SNPsplit.bam
    samtools sort -@ 6 SNPsplit/TW136xTW137_s.genome1.bam -o SNPsplit/TW136xTW137_SNPsplit_g1.bam
    samtools sort -@ 6 SNPsplit/TW136xTW137_s.genome2.bam -o SNPsplit/TW136xTW137_SNPsplit_g2.bam

mkdir featurecounts
ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a nsyl_12.gtf -t 'gene' -g 'ID' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM SNPsplit/TW136xTW137_SNPsplit.bam
featureCounts -T 6 -s 1 -a nsyl_12.gtf -t 'gene' -g 'ID' -o featurecounts/read_counts_g1.tab --readExtension5 500 -R BAM SNPsplit/TW136xTW137_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a nsyl_12.gtf -t 'gene' -g 'ID' -o featurecounts/read_counts_g2.tab --readExtension5 500 -R BAM SNPsplit/TW136xTW137_SNPsplit_g2.bam

conda deactivate