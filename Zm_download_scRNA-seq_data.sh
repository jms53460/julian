#!/bin/bash
#SBATCH --job-name=download_maize_data                  # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=1                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=72:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Zm_download.out         # Location of standard output file
#SBATCH --error=/home/jms53460/Zm_download.err          # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /home/jms53460/Maize_SGT_2022
#curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI10_CKDL190144505-1B-10_HNLFWDSXX_L3_R1.fastq.gz > RPI10_1B-10_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI10_CKDL190144505-1B-10_HNLFWDSXX_L3_R2.fastq.gz > RPI10_1B-10_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI10_TAGCTT_S6_L002_R1.fastq.gz > RPI10_S6_L002_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI10_TAGCTT_S6_L002_R2.fastq.gz > RPI10_S6_L002_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI11_GGCTAC_S5_L002_R1.fastq.gz > RPI11_S5_L002_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI11_GGCTAC_S5_L002_R2.fastq.gz > RPI11_S5_L002_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI1_CKDL190144505-1B-1_HNLFWDSXX_L3_R1.fastq.gz > RPI1_1B-1_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI1_CKDL190144505-1B-1_HNLFWDSXX_L3_R2.fastq.gz > RPI1_1B-1_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI1_CKDL200149798-1a_H7MM2BBXX_L2_R1.fastq.gz > RPI1_1a_L2_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI1_CKDL200149798-1a_H7MM2BBXX_L2_R2.fastq.gz > RPI1_1a_L2_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI2_CKDL190144505-1B-2_HNLFWDSXX_L3_R1.fastq.gz > RPI2_1B-2_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI2_CKDL190144505-1B-2_HNLFWDSXX_L3_R2.fastq.gz > RPI2_1B-2_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI2_CKDL200149799-1a_H7MM2BBXX_L3_R1.fastq.gz > RPI2_1a_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI2_CKDL200149799-1a_H7MM2BBXX_L3_R2.fastq.gz > RPI2_1a_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI3_CKDL190144505-1B-3_HNLFWDSXX_L3_R1.fastq.gz > RPI3_1B-3_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI3_CKDL190144505-1B-3_HNLFWDSXX_L3_R2.fastq.gz > RPI3_1B-3_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI4_CKDL190144505-1B-4_HNLFWDSXX_L3_R1.fastq.gz > RPI4_1B-4_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI4_CKDL190144505-1B-4_HNLFWDSXX_L3_R2.fastq.gz > RPI4_1B-4_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI5_CKDL190144505-1B-5_HNLFWDSXX_L3_R1.fastq.gz > RPI5_1B-5_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI5_CKDL190144505-1B-5_HNLFWDSXX_L3_R2.fastq.gz > RPI5_1B-5_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI6_CKDL190144505-1B-6_HNLFWDSXX_L3_R1.fastq.gz > RPI6_1B-6_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI6_CKDL190144505-1B-6_HNLFWDSXX_L3_R2.fastq.gz > RPI6_1B-6_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI7_CKDL190144505-1B-7_HNLFWDSXX_L3_R1.fastq.gz > RPI7_1B-7_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI7_CKDL190144505-1B-7_HNLFWDSXX_L3_R2.fastq.gz > RPI7_1B-7_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI8_CKDL190144505-1B-8_HNLFWDSXX_L3_R1.fastq.gz > RPI8_1B-8_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI8_CKDL190144505-1B-8_HNLFWDSXX_L3_R2.fastq.gz > RPI8_1B-8_L3_R2.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI9_CKDL190144505-1B-9_HNLFWDSXX_L3_R1.fastq.gz > RPI9_1B-9_L3_R1.fastq.gz
curl -s https://stacks.stanford.edu/file/druid:js115sm3463/RPI9_CKDL190144505-1B-9_HNLFWDSXX_L3_R2.fastq.gz > RPI9_1B-9_L3_R2.fastq.gz
