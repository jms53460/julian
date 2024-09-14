#!/bin/bash
#SBATCH --job-name=Ns_SNPsplit_stringtie                                            # Job name
#SBATCH --partition=batch                                                           # Partition (queue) name
#SBATCH --ntasks=1                                                                  # Single task job
#SBATCH --cpus-per-task=6                                                           # Number of cores per task
#SBATCH --mem=50gb                                                                  # Total memory for job
#SBATCH --time=6:00:00                                                             # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/7_2024_Ns/Ns_SNPsplit_stringtie.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/7_2024_Ns/Ns_SNPsplit_stringtie.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                                # Where to send mail
#SBATCH --mail-type=END,FAIL                                                        # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/7_2024_Ns
#mkdir SNPsplit
#ml SAMtools/1.16.1-GCC-11.3.0
#ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
#for file in "hisat2_out/"*_s.bam
#do
#    file2="${file:11:-6}"
#
#    SNPsplit --conflicting -o SNPsplit --snp_file Ns_SNPs3.tab "$file"
#    samtools sort -@ 6 SNPsplit/"$file2"_s.allele_flagged.bam -o SNPsplit/"$file2"_SNPsplit.bam
    
#done

#for file in "SNPsplit/"*_s.genome1.bam
#do
#    file2="${file:9:-14}"
#    samtools sort -@ 6 "$file" -o SNPsplit/"$file2"_SNPsplit_g1.bam
#done

#for file in "SNPsplit/"*_s.genome2.bam
#do
#    file2="${file:9:-14}"
#    samtools sort -@ 6 "$file" -o SNPsplit/"$file2"_SNPsplit_g2.bam
#done

ml StringTie/2.2.1-GCC-11.3.0
mkdir stringtie_out
for file in "hisat2_out/"*_s.bam
do
	stringtie -p 6 --rf -o "stringtie_out/""${file:11:-6}"".gtf" "$file"
done

# Merge StringTie transcripts
ls -1 "stringtie_out2/"*.gtf | gawk '{print $0}' > mergelist.txt

# Merge GTF files
stringtie --merge -p 6 -o "stringtie_out2/stringtie_merged.gtf" mergelist.txt
rm mergelist.txt
