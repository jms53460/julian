#!/bin/bash
#SBATCH --job-name=Maize_SGT_pipeline                                     # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=400gb                                                       # Total memory for job
#SBATCH --time=48:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_8-10_2025/SGT_Pipeline.out       # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_8-10_2025/SGT_Pipeline.err        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_8-10_2025

mkdir Demultiplexed
ml Mamba/23.11.0-0
source activate /home/jms53460/Fastq-Multx
module load fastp/0.23.4-GCC-13.2.0

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    if [ ! -f "Demultiplexed/""$file2""_1s.fastq.gz" ]; then
	    #Move UMI to header
        fastp -w 6 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed/umi_""$file2""_R1.fastq.gz" -O "Demultiplexed/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI
        
        #Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location
	    fastq-multx -b -B "CELSeq_barcodes.txt" -m 0 "Demultiplexed/umi_""$file2""_R2.fastq.gz" "Demultiplexed/umi_""$file2""_R1.fastq.gz" "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Demultiplexed/""$file2""_%_R2.fastq.gz" "Demultiplexed/""$file2""_%.fastq.gz" "Demultiplexed/""$file2""_%_umi.fastq.gz"

    fi
done
conda deactivate

mkdir SRA_upload
module load fastp/0.23.4-GCC-13.2.0

for file in Demultiplexed/*s.fastq.gz; do
	file2="${file:14:-9}"

    #Trim UMI containing read to only contain the UMI
    fastp -w 6 -B 10 -i "Demultiplexed/""$file2"".fastq.gz" -I "Demultiplexed/""$file2""_umi.fastq.gz" -o "SRA_upload/""$file2"".fastq.gz" -O "SRA_upload/""$file2""_umi.fastq.gz" -A -Q -L -G

done

mkdir hisat2_out
for file in Demultiplexed/*s.fastq.gz; do
	file2="${file:14:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

    #Trim poly A's
	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

ml HISAT2/2.2.1-gompi-2023a
hisat2-build Zm_B73_A188_N-masked_genome.fa Zm_B73_A188_N-masked_genome_index
ml SAMtools/1.21-GCC-13.3.0
for file in hisat2_out/*s.fastq.gz
do
	file2="${file:11:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x Zm_B73_A188_N-masked_genome_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out/""$file2""_s.bam"
	
fi
done


mkdir SNPsplit
module load SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "hisat2_out/"*_s.bam
do
    file2="${file:11:-6}"

    SNPsplit --conflicting -o SNPsplit --snp_file Zm_B73_A188_SNPs.tab "$file"
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


mkdir featurecounts
mkdir bams
mkdir UMIcounts
mkdir UMIcounts_g1
mkdir UMIcounts_g2
ml Mamba/23.11.0-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -t 'gene' -g 'ID' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit.bam
featureCounts -T 6 -s 1 -a Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -t 'gene' -g 'ID' -o featurecounts/read_counts_g1.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -t 'gene' -g 'ID' -o featurecounts/read_counts_g2.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit_g2.bam

conda deactivate


ml SAMtools/1.21-GCC-13.3.0
module load UMI-tools/1.1.4-foss-2023a

for file in "featurecounts/"*SNPsplit.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index -@ 6 "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g1.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g1/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index -@ 6 "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g1/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g2.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g2/${file2}.tsv" ]; then

        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index -@ 6 "bams/$file2"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g2/${file2}.tsv"
    fi
done