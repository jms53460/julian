#!/bin/bash
#SBATCH --job-name=At_features_UMIs2                                          # Job name
#SBATCH --partition=batch                                                     # Partition (queue) name
#SBATCH --ntasks=1                                                            # Single task job
#SBATCH --cpus-per-task=6                                                     # Number of cores per task
#SBATCH --mem=50gb                                                            # Total memory for job
#SBATCH --time=12:00:00                                                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/4_2024_At_Spike_ins/At_features_UMIs2.out  # Location of standard output file
#SBATCH --error=/scratch/jms53460/4_2024_At_Spike_ins/At_features_UMIs2.err   # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                          # Where to send mail
#SBATCH --mail-type=END,FAIL                                                  # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/4_2024_At_Spike_ins
module load SAMtools/1.16.1-GCC-11.3.0
for file in "At_SNPsplit/"*_s.genome1.bam
do
    file2="${file:12:-14}"
    samtools sort -@ 6 "$file" -o At_SNPsplit/"$file2"_SNPsplit_g1.bam
done
for file in "At_SNPsplit/"*_s.genome2.bam
do
    file2="${file:12:-14}"
    samtools sort -@ 6 "$file" -o At_SNPsplit/"$file2"_SNPsplit_g2.bam
done

ml purge_dups/1.2.5-foss-2021b
ml Miniconda3/23.5.2-0
source activate ./subread-env/

featureCounts -T 6 -s 1 -a TAIR10.1_Col_chr.gff -t 'gene' -g 'ID' -o stringtie_out/read_counts_g1.tab --readExtension5 500 -R BAM At_SNPsplit/*_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a stringtie_out/stringtie_merged.gtf -o stringtie_out/read_counts_g1.tab --readExtension5 500 -R BAM At_SNPsplit/*_SNPsplit_g1.bam
featureCounts -T 6 -s 1 -a TAIR10.1_Col_chr.gff -t 'gene' -g 'ID' -o stringtie_out/read_counts_g2.tab --readExtension5 500 -R BAM At_SNPsplit/*_SNPsplit_g2.bam
featureCounts -T 6 -s 1 -a stringtie_out/stringtie_merged.gtf -o stringtie_out/read_counts_g2.tab --readExtension5 500 -R BAM At_SNPsplit/*_SNPsplit_g2.bam

conda deactivate

mkdir UMIcounts_g1
mkdir UMIcounts_g2
for file in "stringtie_out/"*g1.bam*
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

for file in "stringtie_out/"*g2.bam*
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