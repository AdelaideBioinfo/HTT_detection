#!/bin/bash

# Example usage:
# INDIR=results FILE=HTT_condidates.fasta ID=80 sbatch 2_usearchConsensus.sbatch

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=32GB

# go to input dir
cd $INDIR

# run usearch to generate consensus seqs
usearch -cluster_fast $FILE -threads 8 -qmask none -id 0."$ID" -uc "$FILE"_"$ID".uc -consout "$FILE"_"$ID".consensus

