#!/bin/bash
#SBATCH --job-name=MaizeF1_4_2026_hisat_SNPsplit                          # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                       # Total memory for job
#SBATCH --time=12:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/MaizeF1_4_2026/hs_SNP.out              # Location of standard output file
#SBATCH --error=/scratch/jms53460/MaizeF1_4_2026/hs_SNP.err               # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/MaizeF1_4_2026/
mkdir hisat2_out2
ml HISAT2/2.2.1-gompi-2023a
ml SAMtools/1.21-GCC-13.3.0
for file in hisat2_out/*s.fastq.gz
do
	file2="${file:11:-9}"

if [ ! -f "hisat2_out2/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x Zm_B73_A188_N-masked_genome_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out2/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out2/""$file2""_unsorted.bam" -o "hisat2_out2/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out2/""$file2""_s.bam"
	
fi
done

mkdir SNPsplit2
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "hisat2_out2/"*_s.bam
do
    file2="${file:12:-6}"

    SNPsplit --no_sort --single_end --conflicting -o SNPsplit2 --snp_file Zm_B73_A188_SNPs.tab "$file"
    samtools sort -@ 6 SNPsplit2/"$file2"_s.allele_flagged.bam -o SNPsplit2/"$file2"_SNPsplit.bam
    
done

for file in "SNPsplit2/"*_s.genome1.bam
do
    file2="${file:10:-14}"
    samtools sort -@ 6 "$file" -o SNPsplit2/"$file2"_SNPsplit_g1.bam
done

for file in "SNPsplit2/"*_s.genome2.bam
do
    file2="${file:10:-14}"
    samtools sort -@ 6 "$file" -o SNPsplit2/"$file2"_SNPsplit_g2.bam
done