#!/bin/bash

# Example usage:
#
# GENOMEDIR=test_genome \
# GENOME=YarrowiaLipolytica_ASM252v1.fa \
# QUERYDIR=test_query \
# QUERY=L1_nucl_seqs.fasta \
# RESULTSDIR=results \
# sbatch 1b_lastz_and_extract.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=12:00:00 
#SBATCH --mem=32GB    

# Load BEDTools
module load bedtools

echo $GENOME
echo $QUERY

# Run lastz
lastz \
$GENOMEDIR/$GENOME[unmask,multiple] \
$QUERYDIR/$QUERY[unmask,multiple] \
--chain --gapped --coverage=80 --ambiguous=n --ambiguous=iupac \
--format=general-:name2,start2,end2,score,strand2,size2,name1,start1,end1 \
> $GENOME_$QUERY.lastzout

# Rearrange columns to put the concatenated file in BED-like form (for bedtools) 
awk '{print $7 "\t" $8 "\t" $9 "\t" "L1hit" "\t" "1" "\t" $5}' $GENOME_$QUERY.lastzout > bedformat_$GENOME_$QUERY.tmp

# Sort and merge 
bedtools sort -i bedformat_$GENOME_$QUERY.tmp \
| bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct > "$GENOME"_"$QUERY"_lastz.bed

# Extract FASTA from corrected merged BED file
bedtools getfasta -s -fi "$GENOMEDIR/$GENOME" -bed "$GENOME"_"$QUERY"_lastz.bed -fo "$GENOME"_"$QUERY"_lastz.fasta

# Move output bed and fasta to results folder
mv "$GENOME"_"$QUERY"_lastz.fasta $RESULTSDIR
mv "$GENOME"_"$QUERY"_lastz.bed $RESULTSDIR

# Remove temporary files
rm *.tmp
rm *.lastzout
