#!/bin/bash

# Example usage:
#
# INDIR=results/allSpeciesCombined \
# FILE=allSpecies_L1.fasta \
# ID=80 PREFIX=c \
# sbatch 3a_vsearch_cluster_for_nucleotide_seqs.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -c 4
#SBATCH --time=24:00:00 
#SBATCH --mem=32GB  

# go to input dir
cd $INDIR

# make dir for clusters
mkdir -p "$FILE"_"$ID"_clusters

# run vsearch
vsearch -cluster_fast $FILE -threads 4 -qmask none -id 0."$ID" -uc "$FILE"_"$ID".uc -clusters "$FILE"_"$ID"_clusters/"$PREFIX"_