#!/bin/bash
#SBATCH --job-name=Rice_DNA_vcf                                      # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=200gb                                                # Total memory for job
#SBATCH --time=24:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_DNA2/Rice_vcf.out           # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_DNA2/Rice_vcf.err            # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_DNA2
module load SAMtools/1.21-GCC-13.3.0
for file in "Raw_Data/"*_R1_001.fastq.gz
do
	file2="${file:9:-16}"
samtools index -@ 6 hisat2_out/"$file2"_s.bam

module load BCFtools/1.21-GCC-13.3.0
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f Nipponbare_12.fna hisat2_out/"$file2"_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > "$file2".vcf.gz
bcftools index "$file2".vcf.gz
done

gzip -d "$file2".vcf.gz
bcftools view -i 'TYPE="snp"' "$file2".vcf > "$file2"_snps.vcf
bcftools view -i 'GT="het"' "$file2"_snps.vcf > "$file2"_het_snps.vcf

###Now I will convert the vcf to a tsv file and edit it to be a SNP file ready for SNPsplit
ml Mamba/23.11.0-0
source activate /home/jms53460/vcf2tsvpy
vcf2tsvpy --input_vcf "$file2"_het_snps.vcf --out_tsv "$file2"_vcf_table.tsv 
conda deactivate

###Selected columns from the vcf_table
awk '{print $3,$1,$2,$6,$4,$5}' "$file2"_vcf_table.tsv OFS="\t" > "$file2"_snps.tsv