#!/bin/bash
#SBATCH --job-name=SNPsplit                                                      # Job name
#SBATCH --partition=batch                                                           # Partition (queue) name
#SBATCH --ntasks=1                                                                  # Single task job
#SBATCH --cpus-per-task=6                                                           # Number of cores per task
#SBATCH --mem=70gb                                                                  # Total memory for job
#SBATCH --time=12:00:00                                                              # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/12_2024_Sl/SNPsplit.out                       # Location of standard output file
#SBATCH --error=/scratch/jms53460/12_2024_Sl/SNPsplit.err                        # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                                # Where to send mail
#SBATCH --mail-type=END,FAIL                                                        # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/12_2024_Sl
mkdir SNPsplit2
ml SAMtools/1.16.1-GCC-11.3.0
ml SNPsplit/0.6.0-GCC-11.3.0-Perl-5.34.1
for file in "hisat2_out2/"*_s.bam
do
    file2="${file:12:-6}"

    SNPsplit --conflicting -o SNPsplit2 --snp_file Sl_SNPs.tab "$file"
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

mkdir featurecounts2
mkdir bams2
mkdir UMIcounts2
mkdir UMIcounts2_g1
mkdir UMIcounts2_g2
ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a Sl_14.gff -t 'gene' -g 'ID' -o featurecounts2/read_counts.tab --readExtension5 500 -R BAM SNPsplit2/*_SNPsplit.bam
featureCounts -T 6 -s 1 -a Sl_14.gff -t 'gene' -g 'ID' -o featurecounts2/read_counts_g1.tab --readExtension5 500 -R BAM SNPsplit2/*_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a Sl_14.gff -t 'gene' -g 'ID' -o featurecounts2/read_counts_g2.tab --readExtension5 500 -R BAM SNPsplit2/*_SNPsplit_g2.bam

conda deactivate

for file in "featurecounts2/"*SNPsplit.bam*
do
    file2="${file:15:-22}"
    if [ ! -f "UMIcounts2/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams2/$file2"
        samtools index "bams2/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams2/$file2" -S "UMIcounts2/${file2}.tsv"
    fi
done

for file in "featurecounts2/"*g1.bam*
do
    file2="${file:15:-22}"
    if [ ! -f "UMIcounts2_g1/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams2/$file2"
        samtools index "bams2/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams2/$file2" -S "UMIcounts2_g1/${file2}.tsv"
    fi
done

for file in "featurecounts2/"*g2.bam*
do
    file2="${file:15:-22}"
    if [ ! -f "UMIcounts2_g2/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams2/$file2"
        samtools index "bams2/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams2/$file2" -S "UMIcounts2_g2/${file2}.tsv"
    fi
done
