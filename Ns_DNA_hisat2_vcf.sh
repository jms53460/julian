#!/bin/bash
#SBATCH --job-name=Ns_hisat2_vcf                                      # Job name
#SBATCH --partition=batch                                             # Partition (queue) name
#SBATCH --ntasks=1                                                    # Single task job
#SBATCH --cpus-per-task=6                                             # Number of cores per task
#SBATCH --mem=50gb                                                    # Total memory for job
#SBATCH --time=12:00:00                                               # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_h2_vcf.out       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_h2_vcf.err        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                  # Where to send mail
#SBATCH --mail-type=END,FAIL                                          # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
mkdir hisat2_out2
cp hisat2_out/*fastq.gz hisat2_out2
for file in "Raw_Data/"*_R1_001.fastq.gz
do
	file2="${file:9:-16}"

if [ ! -f "hisat2_out2/""$file2"".bam" ]; then

    ml HISAT2/3n-20201216-gompi-2022a
    hisat2-build TW136_consensus.fa TW136_consensus_index
	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 6 --dta -x TW136_consensus_index -1 "hisat2_out2/""$file2""_R1_001.fastq.gz" -2 "hisat2_out2/""$file2""_R2_001.fastq.gz" | samtools view -bS -> "hisat2_out2/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out2/""$file2""_unsorted.bam" -o "hisat2_out2/""$file2""_s.bam"
	
fi
done

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f TW136_consensus.fa hisat2_out2/TW137_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW137_consensus.vcf.gz
bcftools index TW137_consensus.vcf.gz
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f TW136_consensus.fa hisat2_out2/TW136_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW136_consensus.vcf.gz
bcftools index TW136_consensus.vcf.gz