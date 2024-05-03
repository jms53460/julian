#!/bin/bash
#SBATCH --job-name=At_SNPsplit_stringtie                                            # Job name
#SBATCH --partition=batch                                                           # Partition (queue) name
#SBATCH --ntasks=1                                                                  # Single task job
#SBATCH --cpus-per-task=6                                                           # Number of cores per task
#SBATCH --mem=50gb                                                                  # Total memory for job
#SBATCH --time=12:00:00                                                             # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/4_2024_At_Spike_ins/At_SNPsplit_stringtie.out    # Location of standard output file
#SBATCH --error=/scratch/jms53460/4_2024_At_Spike_ins/At_SNPsplit_stringtie.err     # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                                # Where to send mail
#SBATCH --mail-type=END,FAIL                                                        # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/4_2024_At_Spike_ins
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "hisat2_out/"*_s.bam
do
    file2="${file:11:-6}"

    SNPsplit --conflicting -o At_SNPsplit --snp_file Ler_SNPs.tab "$file"
    samtools sort -@ 6 At_SNPsplit/"$file2"_s.allele_flagged.bam -o At_SNPsplit/"$file2"_SNPsplit.bam
    
done

ml StringTie/2.2.1-GCC-11.3.0
for file in "hisat2_out/"*_s.bam
do
    file2="${file:11:-6}"

    stringtie At_SNPsplit/"$file2"_SNPsplit.bam -p 6 -G TAIR10.1_Col_chr.gff --rf -o stringtie_out/"$file2".gtf
done    
    
#Merge stringtie transcripts
ls -1 "stringtie_out/"*.gtf | gawk '{print $0}' > mergelist.txt
# Merge GTF files
stringtie --merge -p 6 -G TAIR10.1_Col_chr.gff -o "stringtie_out/stringtie_merged.gtf" mergelist.txt
