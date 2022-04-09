#!/bin/bash

# Want to keep only the nucleotide sequences that CENSOR recognises as your repeat element of interest (e.g. L1s)
# And not some other repeat (e.g. CR1 or hAT)

# Example usage
# 
# SPECIES=Yarrowia.lipolytica \
# FILE=Yarrowia.lipolytica_L1_combined.fasta \
# GENOME=test_genome/YarrowiaLipolytica_ASM252v1.fa \
# ELEMENT=L1 
# QUERY=test_query/known_L1_elements_from_repbase.txt \
# CENSORDIR=results/censored \
# sbatch 2b_check_censor_output.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --time=12:00:00 
#SBATCH --mem=32GB

# Load the necessary modules
module load bedtools

# First, extract map file lines in the right orientation (d, not c)
# Print the header (which contains coords and strand) and column 4 (hit name)
# Put into BED-like format
cat $CENSORDIR/"$FILE".map \
| awk '{if ($7=="d") print $0 }' \
| awk '{print $1 "\t" $4}' \
| sed 's/:/\t/g' \
| sed 's/(-)/.rev/g' \
| sed 's/(+)/.fwd/g' \
| awk '{gsub(/-/,"\t",$2); print}' \
| sed 's/.rev/\t-/g' \
| sed 's/.fwd/\t+/g' \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' \
> "$FILE".rearranged.tmp

# Sort, merge and rearrange columns
bedtools sort -i "$FILE".rearranged.tmp \
| bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct \
| awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $6}' \
> "$FILE".merged.tmp

# Extract sequences which correctly identified as L1s
awk 'NR==FNR { a[$1] = $1; next} { for (k in a) if ($1 ~ a[k]) { print $0; break } }' "$QUERY" "$FILE".merged.tmp \
| awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5}' \
> "$SPECIES"_"$ELEMENT"_censored.bed

# Extract FASTA from L1-only merged BED file
bedtools getfasta -s -fi "$GENOME" -bed "$SPECIES"_"$ELEMENT"_censored.bed -fo "$SPECIES"_"$ELEMENT"_censored.tmp

# Sort sequences by length
usearch -sortbylength "$SPECIES"_"$ELEMENT"_censored.tmp -fastaout "$SPECIES"_"$ELEMENT"_censored.fasta

# Move files to censor dir
mv "$SPECIES"_"$ELEMENT"_censored.bed $CENSORDIR
mv "$SPECIES"_"$ELEMENT"_censored.fasta $CENSORDIR

# Remove unnecessary files
rm *.tmp
