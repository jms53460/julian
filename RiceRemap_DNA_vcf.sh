#!/bin/bash
#SBATCH --job-name=RiceRemap_DNA_vcf                                      # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=200gb                                                # Total memory for job
#SBATCH --time=24:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/RiceRemap/Rice_vcf.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/RiceRemap/Rice_vcf.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/RiceRemap
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
for file in "Raw_Data_DNA/"*_R1_001.fastq.gz
do
	file2="${file:13:-16}"
samtools index -@ 6 hisat2_out_DNA/"$file2"_s.bam

bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f Rice_r7.fa hisat2_out_DNA/"$file2"_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > "$file2".vcf.gz
bcftools index "$file2".vcf.gz
done