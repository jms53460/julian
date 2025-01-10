#!/bin/bash
#SBATCH --job-name=Ns_DNA_vcf2                                      # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=100gb                                                 # Total memory for job
#SBATCH --time=24:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_vcf2.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_vcf2.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
module load BCFtools/1.15.1-GCC-11.3.0
# normalize indels
bcftools norm -f Ns_genome.fna TW136.vcf.gz -Ob -o TW136.norm.bcf 

# filter adjacent indels within 5bp
bcftools filter --IndelGap 5 TW136.norm.bcf -Ob -o TW136.norm.flt-indels.bcf

# apply variants to create consensus sequence
cat Ns_genome.fna | bcftools consensus TW136.vcf.gz > TW136_consensus.fa

# output IUPAC ambiguity codes based on REF+ALT columns (regardless of genotype)
cat Ns_genome.fna | bcftools consensus --iupac-codes TW136.vcf.gz > TW136_consensus.fa

# output IUPAC ambiguity codes based on sample genotypes
cat Ns_genome.fna | bcftools consensus --haplotype I TW136.vcf.gz > TW136_consensus.fa

bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f TW136_consensus.fa hisat2_out/TW137_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW137_consensus.vcf.gz
bcftools index TW137_consensus.vcf.gz
