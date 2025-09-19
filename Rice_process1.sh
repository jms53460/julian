#!/bin/bash
#SBATCH --job-name=Rice_process1                                          # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=100gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_8_2025/NRE1-1/Hs2.out             # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_8_2025/NRE1-1/Hs2.err              # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_8_2025/NRE1-1

module load fastp/0.23.4-GCC-13.2.0
mkdir hisat2_out
for file in Demultiplexed/*s.fastq.gz; do
	file2="${file:14:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	fastp -w 6 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

ml HISAT2/2.2.1-gompi-2023a
hisat2-build NRE1-1_N-masked_12.fa NRE1-1_N-masked_12_index
ml SAMtools/1.21-GCC-13.3.0
for file in hisat2_out/*s.fastq.gz
do
	file2="${file:11:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

	hisat2 -p 6 --dta -x NRE1-1_N-masked_12_index -U "hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "hisat2_out/""$file2""_unsorted.bam"
	samtools sort -@ 6 "hisat2_out/""$file2""_unsorted.bam" -o "hisat2_out/""$file2""_s.bam"
    samtools index -@ 6 "hisat2_out/""$file2""_s.bam"
	
fi
done


mkdir SNPsplit
ml Mamba/23.11.0-0
conda activate /home/jms53460/SNPsplit_env
for file in "hisat2_out/"*_s.bam
do
    file2="${file:11:-6}"

    SNPsplit --conflicting -o SNPsplit --snp_file NRE1-1_SNPs.tab "$file"
    samtools sort -@ 6 SNPsplit/"$file2"_s.allele_flagged.bam -o SNPsplit/"$file2"_SNPsplit.bam
    
done
conda deactivate

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

source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a Nipponbare_12.gtf -t 'gene' -g 'gene_id' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit.bam
featureCounts -T 6 -s 1 -a Nipponbare_12.gtf -t 'gene' -g 'gene_id' -o featurecounts/read_counts_g1.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a Nipponbare_12.gtf -t 'gene' -g 'gene_id' -o featurecounts/read_counts_g2.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit_g2.bam

conda deactivate



for file in "featurecounts/"*SNPsplit.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        ml SAMtools/1.21-GCC-13.3.0
        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        module load UMI-tools/1.1.4-foss-2023a
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g1.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g1/${file2}.tsv" ]; then

        ml SAMtools/1.21-GCC-13.3.0
        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        module load UMI-tools/1.1.4-foss-2023a
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g1/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g2.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g2/${file2}.tsv" ]; then

        ml SAMtools/1.21-GCC-13.3.0
        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        module load UMI-tools/1.1.4-foss-2023a
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g2/${file2}.tsv"
    fi
done