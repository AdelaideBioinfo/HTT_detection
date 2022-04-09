#!/bin/bash

# Example usage:
#
# INDIR=results/allSpeciesCombined \
# FILE=allSpecies_L1_ORFp.fasta \
# ID=90 PREFIX=r \
# sbatch 3b_usearch_cluster_for_aa_seqs.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -c 4
#SBATCH --time=24:00:00 
#SBATCH --mem=32GB  

# go to input dir
cd $INDIR

# make dir for clusters
mkdir -p "$FILE"_"$ID"_clusters

# run usearch
# important note: the open source version (vsearch) does not support protein sequences
# so have to use usearch for any amino acid sequence clustering
usearch -cluster_fast $FILE -threads 4 -qmask none -id 0."$ID" -uc "$FILE"_"$ID".uc -clusters "$FILE"_"$ID"_clusters/"$PREFIX"_
