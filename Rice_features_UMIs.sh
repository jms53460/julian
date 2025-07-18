#!/bin/bash
#SBATCH --job-name=Features_UMIs                                              # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Rice_7_2025/Features_UMIs.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/Rice_7_2025/Features_UMIs.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Rice_7_2025
mkdir featurecounts
mkdir bams
mkdir UMIcounts
mkdir UMIcounts_g1
mkdir UMIcounts_g2
ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate /home/jms53460/subread-env

featureCounts -T 6 -s 1 -a Nipponbare_12.gtf -t 'gene' -g 'gene_id' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit.bam
featureCounts -T 6 -s 1 -a Nipponbare_12.gtf -t 'gene' -g 'gene_id' -o featurecounts/read_counts_g1.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a Nipponbare_12.gtf -t 'gene' -g 'gene_id' -o featurecounts/read_counts_g2.tab --readExtension5 500 -R BAM SNPsplit/*_SNPsplit_g2.bam

conda deactivate

for file in "featurecounts/"*SNPsplit.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g1.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g1/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g1/${file2}.tsv"
    fi
done

for file in "featurecounts/"*g2.bam*
do
    file2="${file:14:-22}"
    if [ ! -f "UMIcounts_g2/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 6 "$file" -o "bams/$file2"
        samtools index "bams/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "bams/$file2" -S "UMIcounts_g2/${file2}.tsv"
    fi
done