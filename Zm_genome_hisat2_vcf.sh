#!/bin/bash
#SBATCH --job-name=Zm_genome_hisat2_vcf                                    # Job name
#SBATCH --partition=batch                                                  # Partition (queue) name
#SBATCH --ntasks=1                                                         # Single task job
#SBATCH --cpus-per-task=6                                                  # Number of cores per task
#SBATCH --mem=50gb                                                         # Total memory for job
#SBATCH --time=12:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_genome_hisat2_vcf.out # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_genome_hisat2_vcf.err  # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                       # Where to send mail
#SBATCH --mail-type=END,FAIL                                               # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022

ml wgsim/20111017-GCC-11.3.0
wgsim -1 10000 -2 10000 -d 20000 A188_genome.fa A188_sim.read1.fq A188_sim.read2.fq 
gzip A188_sim.read1.fq
gzip A188_sim.read2.fq

ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0
hisat2 -p 6 --dta -x B73_v5_hisat2_index -1 A188_sim.read1.fq.gz -2 A188_sim.read2.fq.gz | samtools view -bS -> hisat2_out/A188_B73_unsorted.bam
samtools sort -@ 6 hisat2_out/A188_B73_unsorted.bam -o A188_B73_s.bam
samtools index -@ 6 A188_B73_s.bam

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou -d 1000000 --threads 6 --min-MQ 60 -f B73_v5_genome.fa A188_B73_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > Zm_A188_B73_vcf.gz
bcftools index Zm_A188_B73_vcf.gz