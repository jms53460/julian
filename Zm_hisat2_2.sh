#!/bin/bash
#SBATCH --job-name=Zm_hisat2_2                                            # Job name
#SBATCH --partition=batch                                                 # Partition (queue) name
#SBATCH --ntasks=1                                                        # Single task job
#SBATCH --cpus-per-task=6                                                 # Number of cores per task
#SBATCH --mem=50gb                                                        # Total memory for job
#SBATCH --time=24:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=/scratch/jms53460/Maize_SGT_2022/Zm_hs2_2.out            # Location of standard output file
#SBATCH --error=/scratch/jms53460/Maize_SGT_2022/Zm_hs2_2.err             # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                                      # Where to send mail
#SBATCH --mail-type=END,FAIL                                              # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/jms53460/Maize_SGT_2022
ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0
hisat2 -p 6 --dta -x /scratch/jms53460/Maize_SGT_2022/Zm_N-masked_hisat2_index -U hisat2_out2/RPI1_1a_L2_10s.fastq.gz,hisat2_out2/RPI1_1a_L2_11s.fastq.gz,hisat2_out2/RPI1_1a_L2_12s.fastq.gz,hisat2_out2/RPI1_1a_L2_13s.fastq.gz,hisat2_out2/RPI1_1a_L2_14s.fastq.gz,hisat2_out2/RPI1_1a_L2_15s.fastq.gz,hisat2_out2/RPI1_1a_L2_16s.fastq.gz,hisat2_out2/RPI1_1a_L2_17s.fastq.gz,hisat2_out2/RPI1_1a_L2_18s.fastq.gz,hisat2_out2/RPI1_1a_L2_19s.fastq.gz,hisat2_out2/RPI1_1a_L2_1s.fastq.gz,hisat2_out2/RPI1_1a_L2_20s.fastq.gz,hisat2_out2/RPI1_1a_L2_21s.fastq.gz,hisat2_out2/RPI1_1a_L2_22s.fastq.gz,hisat2_out2/RPI1_1a_L2_23s.fastq.gz,hisat2_out2/RPI1_1a_L2_24s.fastq.gz,hisat2_out2/RPI1_1a_L2_25s.fastq.gz,hisat2_out2/RPI1_1a_L2_26s.fastq.gz,hisat2_out2/RPI1_1a_L2_27s.fastq.gz,hisat2_out2/RPI1_1a_L2_28s.fastq.gz,hisat2_out2/RPI1_1a_L2_29s.fastq.gz,hisat2_out2/RPI1_1a_L2_2s.fastq.gz,hisat2_out2/RPI1_1a_L2_30s.fastq.gz,hisat2_out2/RPI1_1a_L2_31s.fastq.gz,hisat2_out2/RPI1_1a_L2_32s.fastq.gz,hisat2_out2/RPI1_1a_L2_33s.fastq.gz,hisat2_out2/RPI1_1a_L2_34s.fastq.gz,hisat2_out2/RPI1_1a_L2_35s.fastq.gz,hisat2_out2/RPI1_1a_L2_36s.fastq.gz,hisat2_out2/RPI1_1a_L2_37s.fastq.gz,hisat2_out2/RPI1_1a_L2_38s.fastq.gz,hisat2_out2/RPI1_1a_L2_39s.fastq.gz,hisat2_out2/RPI1_1a_L2_3s.fastq.gz,hisat2_out2/RPI1_1a_L2_40s.fastq.gz,hisat2_out2/RPI1_1a_L2_41s.fastq.gz,hisat2_out2/RPI1_1a_L2_42s.fastq.gz,hisat2_out2/RPI1_1a_L2_43s.fastq.gz,hisat2_out2/RPI1_1a_L2_44s.fastq.gz,hisat2_out2/RPI1_1a_L2_45s.fastq.gz,hisat2_out2/RPI1_1a_L2_46s.fastq.gz,hisat2_out2/RPI1_1a_L2_47s.fastq.gz,hisat2_out2/RPI1_1a_L2_48s.fastq.gz,hisat2_out2/RPI1_1a_L2_49s.fastq.gz,hisat2_out2/RPI1_1a_L2_4s.fastq.gz,hisat2_out2/RPI1_1a_L2_50s.fastq.gz,hisat2_out2/RPI1_1a_L2_51s.fastq.gz,hisat2_out2/RPI1_1a_L2_52s.fastq.gz,hisat2_out2/RPI1_1a_L2_53s.fastq.gz,hisat2_out2/RPI1_1a_L2_54s.fastq.gz,hisat2_out2/RPI1_1a_L2_55s.fastq.gz,hisat2_out2/RPI1_1a_L2_56s.fastq.gz,hisat2_out2/RPI1_1a_L2_57s.fastq.gz,hisat2_out2/RPI1_1a_L2_58s.fastq.gz,hisat2_out2/RPI1_1a_L2_59s.fastq.gz,hisat2_out2/RPI1_1a_L2_5s.fastq.gz,hisat2_out2/RPI1_1a_L2_60s.fastq.gz,hisat2_out2/RPI1_1a_L2_61s.fastq.gz,hisat2_out2/RPI1_1a_L2_62s.fastq.gz,hisat2_out2/RPI1_1a_L2_63s.fastq.gz,hisat2_out2/RPI1_1a_L2_64s.fastq.gz,hisat2_out2/RPI1_1a_L2_65s.fastq.gz,hisat2_out2/RPI1_1a_L2_66s.fastq.gz,hisat2_out2/RPI1_1a_L2_67s.fastq.gz,hisat2_out2/RPI1_1a_L2_68s.fastq.gz,hisat2_out2/RPI1_1a_L2_69s.fastq.gz,hisat2_out2/RPI1_1a_L2_6s.fastq.gz,hisat2_out2/RPI1_1a_L2_70s.fastq.gz,hisat2_out2/RPI1_1a_L2_71s.fastq.gz,hisat2_out2/RPI1_1a_L2_72s.fastq.gz,hisat2_out2/RPI1_1a_L2_73s.fastq.gz,hisat2_out2/RPI1_1a_L2_74s.fastq.gz,hisat2_out2/RPI1_1a_L2_75s.fastq.gz,hisat2_out2/RPI1_1a_L2_76s.fastq.gz,hisat2_out2/RPI1_1a_L2_77s.fastq.gz,hisat2_out2/RPI1_1a_L2_78s.fastq.gz,hisat2_out2/RPI1_1a_L2_79s.fastq.gz,hisat2_out2/RPI1_1a_L2_7s.fastq.gz,hisat2_out2/RPI1_1a_L2_80s.fastq.gz,hisat2_out2/RPI1_1a_L2_81s.fastq.gz,hisat2_out2/RPI1_1a_L2_82s.fastq.gz,hisat2_out2/RPI1_1a_L2_83s.fastq.gz,hisat2_out2/RPI1_1a_L2_84s.fastq.gz,hisat2_out2/RPI1_1a_L2_85s.fastq.gz,hisat2_out2/RPI1_1a_L2_86s.fastq.gz,hisat2_out2/RPI1_1a_L2_87s.fastq.gz,hisat2_out2/RPI1_1a_L2_88s.fastq.gz,hisat2_out2/RPI1_1a_L2_89s.fastq.gz,hisat2_out2/RPI1_1a_L2_8s.fastq.gz,hisat2_out2/RPI1_1a_L2_90s.fastq.gz,hisat2_out2/RPI1_1a_L2_91s.fastq.gz,hisat2_out2/RPI1_1a_L2_92s.fastq.gz,hisat2_out2/RPI1_1a_L2_93s.fastq.gz,hisat2_out2/RPI1_1a_L2_94s.fastq.gz,hisat2_out2/RPI1_1a_L2_95s.fastq.gz,hisat2_out2/RPI1_1a_L2_96s.fastq.gz,hisat2_out2/RPI1_1a_L2_9s.fastq.gz,hisat2_out2/RPI1_1a_L2_unmatched.fastq.gz | samtools view -bS -> hisat2_out2/merged_unsorted.bam
samtools sort -@ 6 hisat2_out2/merged_unsorted.bam -o hisat2_out2/merged_s.bam
samtools index -@ 6 hisat2_out2/merged_s.bam