#!/bin/bash
#SBATCH --job-name=Ns_nmask                                               # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=12:00:00                                                    # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Tobacco_DNAseq/Ns_nmask.out                # Location of standard output file
#SBATCH --error=/scratch/jms53460/Tobacco_DNAseq/Ns_nmask.err                 # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Tobacco_DNAseq
ml BEDTools/2.30.0-GCC-12.2.0
bedtools maskfasta -fi TW136_consensus.fa -fo TW136_consensus_N-masked.fa -bed TW137_hom_SNPs.vcf -fullHeader

module load BCFtools/1.15.1-GCC-11.3.0
bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f TW136_consensus.fa hisat2_out2/TW136_s.bam | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > TW136_consensus.vcf.gz
bcftools index TW136_consensus.vcf.gz