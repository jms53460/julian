#!/bin/bash
#SBATCH --job-name=RiceRemap_process3                                          # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=200gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/RiceRemap/NRE1-3/Hs2.out             # Location of standard output file
#SBATCH --error=/scratch/jms53460/RiceRemap/NRE1-3/Hs2.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/RiceRemap/NRE1-3

module load fastp/0.23.4-GCC-13.2.0
mkdir hisat2_out
for file in Demultiplexed/*s.fastq.gz; do
	file2="${file:14:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

ml HISAT2/2.2.1-gompi-2023a
hisat2-build NRE1-3_N-masked_Remap.fa NRE1-3_N-masked_Remap_index
ml SAMtools/1.21-GCC-13.3.0
for file in hisat2_out/*s.fastq.gz
do
	file2="${file:11:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x NRE1-3_N-masked_Remap_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out/""$file2""_s.bam"
	
fi
done

mkdir SNPsplit
module load SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "hisat2_out/"*_s.bam
do
    file2="${file:11:-6}"

    SNPsplit --conflicting -o SNPsplit --snp_file NRE1-3_SNPs_Remap.tab "$file"
    samtools sort -@ 6 SNPsplit/"$file2"_s.allele_flagged.bam -o SNPsplit/"$file2"_SNPsplit.bam
    
done

for file in "SNPsplit/"*_s.genome1.bam
do
    file2="${file:9:-14}"
    samtools sort -@ 6 "$file" -o SNPsplit/"$file2"_SNPsplit_g1.bam
done

for file in "SNPsplit/"*_s.genome2.bam
do
    file2="${file:9:-14}"
    samtools sort -@ 6 "$file" -o SNPsplit/"$file2"_SNPsplit_g2.bam
done