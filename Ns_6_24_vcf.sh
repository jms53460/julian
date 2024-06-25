#!/bin/bash
#SBATCH --job-name=Ns_6_24_vcf                                     # Job name
#SBATCH --partition=batch                                          # Partition (queue) name
#SBATCH --ntasks=1                                                 # Single task job
#SBATCH --cpus-per-task=6                                          # Number of cores per task
#SBATCH --mem=50gb                                                 # Total memory for job
#SBATCH --time=12:00:00                                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/June2024Seq/Ns_vcf.out          # Location of standard output file
#SBATCH --error=/scratch/jms53460/June2024Seq/Ns_vcf.err           # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                               # Where to send mail
#SBATCH --mail-type=END,FAIL                                       # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/June2024Seq
module load SAMtools/1.16.1-GCC-11.3.0

for file in "hisat2_out/"*_s.bam
do
	file2="${file:15:-19}"

    samtools index -@ 6 "$file"

    module load BCFtools/1.15.1-GCC-11.3.0
    bcftools mpileup -Ou --threads 6 -d 1000000 --min-MQ 60 -f Ns_genome.fna "$file" | bcftools call -Ou -m -v --threads 6 | bcftools filter -Oz -e 'QUAL<40 || DP<10' > "$file2".vcf.gz
    bcftools index "$file2".vcf.gz
done