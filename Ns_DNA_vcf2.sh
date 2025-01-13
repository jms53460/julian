#!/bin/bash
#SBATCH --job-name=Ns_DNA_vcf2                                      # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=100gb                                                 # Total memory for job
#SBATCH --time=12:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_vcf2.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_vcf2.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
module load BCFtools/1.15.1-GCC-11.3.0

# apply variants to create consensus sequence
cat Ns_genome.fna | bcftools consensus TW136_hom_SNPs1.vcf.gz > TW136_consensus.fa

# output IUPAC ambiguity codes based on REF+ALT columns (regardless of genotype)
cat Ns_genome.fna | bcftools consensus --iupac-codes TW136_hom_SNPs1.vcf.gz > TW136_consensus.fa

# output IUPAC ambiguity codes based on sample genotypes
cat Ns_genome.fna | bcftools consensus --haplotype I TW136_hom_SNPs1.vcf.gz > TW136_consensus.fa

