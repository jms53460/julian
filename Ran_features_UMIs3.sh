#!/bin/bash
#SBATCH --job-name=Ran_features_UMIs                                      # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=1                                                 # Number of cores per task
#SBATCH --mem=600gb                                                       # Total memory for job
#SBATCH --time=48:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Ran/Ran_features_UMIs3.out              # Location of standard output file
#SBATCH --error=/scratch/jms53460/Ran/Ran_features_UMIs3.err               # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Ran

#mkdir featurecounts
#mkdir bams
#mkdir UMIcounts
#ml Mamba/23.11.0-0
#source activate /home/jms53460/subread-env

#featureCounts -T 6 -s 1 -a JEC21_56.gff -t 'protein_coding_gene' -g 'ID' -o featurecounts/read_counts.tab --readExtension5 500 -R BAM Mapped/*.bam

#conda deactivate


#ml SAMtools/1.21-GCC-13.3.0
#for file in "featurecounts2/"*.bam
#do
#    file2="${file:15:-22}"

#        samtools sort -@ 6 "$file" -o "bams2/$file2"
#        samtools index "bams2/$file2"
#done


module load UMI-tools/1.1.4-foss-2023a
for file in "bams5/"*.bam
do
    file2="${file:6:-25}"

        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "$file" -S "UMIcounts3/"Ran""${file2}.tsv"
done