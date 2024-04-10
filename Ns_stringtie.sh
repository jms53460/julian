#!/bin/bash
#SBATCH --job-name=Ns_stringtie                         # Job name
#SBATCH --partition=batch                               # Partition (queue) name
#SBATCH --ntasks=1                                      # Single task job
#SBATCH --cpus-per-task=6                               # Number of cores per task
#SBATCH --mem=50gb                                      # Total memory for job
#SBATCH --time=12:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=/home/jms53460/Ns_stringtie.out        # Location of standard output file
#SBATCH --error=/home/jms53460/Ns_stringtie.err         # Location of error log file
#SBATCH --mail-user=jms53460@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

for file in "hisat2_out/"*_s.bam
do
	ml StringTie/2.2.1-GCC-11.3.0
	stringtie -p 6 --rf -o "stringtie_out2/""${file:11:-4}"".gtf" "$file"
done

# Merge StringTie transcripts
ls -1 "stringtie_out2/"*.gtf | gawk '{print $0}' > mergelist.txt

# Merge GTF files
stringtie --merge -p 6 -o "stringtie_out2/stringtie_merged.gtf" mergelist.txt
rm mergelist.txt