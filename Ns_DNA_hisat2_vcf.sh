#!/bin/bash
#SBATCH --job-name=Ns_hisat2_vcf                                      # Job name
#SBATCH --partition=batch                                             # Partition (queue) name
#SBATCH --ntasks=1                                                    # Single task job
#SBATCH --cpus-per-task=6                                             # Number of cores per task
#SBATCH --mem=50gb                                                    # Total memory for job
#SBATCH --time=24:00:00                                               # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_h2_vcf.out       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_h2_vcf.err        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                  # Where to send mail
#SBATCH --mail-type=END,FAIL                                          # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
mkdir hisat2_out4
cp hisat2_out/*fastq.gz hisat2_out4

    ml HISAT2/3n-20201216-gompi-2022a
    hisat2-build TW136_consensus2.fa TW136_consensus2_index
	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 6 --dta -x TW136_consensus2_index -1 "hisat2_out4/TW136_R1_001.fastq.gz" -2 "hisat2_out4/TW136_R2_001.fastq.gz" | samtools view -bS -> "hisat2_out4/TW136_unsorted.bam"
	samtools sort -@ 6 "hisat2_out4/TW136_unsorted.bam" -o "hisat2_out4/TW136_s.bam"
	hisat2 -p 6 --dta -x TW136_consensus2_index -1 "hisat2_out4/TW137_R1_001.fastq.gz" -2 "hisat2_out4/TW137_R2_001.fastq.gz" | samtools view -bS -> "hisat2_out4/TW137_unsorted.bam"
	samtools sort -@ 6 "hisat2_out4/TW137_unsorted.bam" -o "hisat2_out4/TW137_s.bam"


module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f TW136_consensus2.fa hisat2_out4/TW136_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW136_consensus2.vcf.gz
bcftools index TW136_consensus2.vcf.gz
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f TW136_consensus2.fa hisat2_out4/TW137_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW137_consensus2.vcf.gz
bcftools index TW137_consensus2.vcf.gz
