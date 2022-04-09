#!/bin/bash

# Example usage:
#
# SPECIES=YarrowiaLipolytica \
# ELEMENT=L1 \
# LASTZFILE=YarrowiaLipolytica_ASM252v1.fa_L1_nucl_seqs.fasta_lastz.bed \
# TBLASTNFILE=YarrowiaLipolytica_ASM252v1.fa_L1_ORFp.fasta_merged.bed \
# GENOME=YarrowiaLipolytica_ASM252v1.fa \
# RESULTSDIR=results \
# sbatch 1c_combine_hits.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=12:00:00 
#SBATCH --mem=32GB    

# Load necessary modules
module load bedtools

# Combine bed files
cat $RESULTSDIR/$TBLASTNFILE $RESULTSDIR/$LASTZFILE > "$SPECIES"_"$ELEMENT"_together.tmp

# Sort and merge hits
bedtools sort -i "$SPECIES"_"$ELEMENT"_together.tmp \
| bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct > $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.bed

# Extract nucleotide sequences from bed file
bedtools getfasta -s -fi test_genome/$GENOME -bed $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.bed -fo $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.fasta

# Remove tmp files
rm *tmp
