#!/bin/bash

# Example usage:
#
# INDIR=results FILE=YarrowiaLipolytica_L1_combined.fasta OUTDIR=results/censored sbatch 2a_censor_sequences.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -c 8
#SBATCH --time=24:00:00 
#SBATCH --mem=32GB

# Load the necessary modules
module load bioperl
module load censor
module load wu-blast
module load bedtools

# Run CENSOR on the extracted nucleotide seqs, using the Repbase library
# Note: download the latest Repbase (http://www.girinst.org/repbase) or Dfam (https://www.dfam.org/releases) repeat library
censor -bprm cpus=8 -lib RepBase21.03_all_seqs.ref "$INDIR/$FILE"

# Move .map file to results directory
mv "$FILE".map $OUTDIR

# Remove or move temporary files
rm "$FILE".*
rm censor.*.log
