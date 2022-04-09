#!/bin/bash

# Example usage:
#
# DIR=test_genome \
# DATABASE=YarrowiaLipolytica_ASM252v1.fa \
# QUERY=L1_ORFp.fasta \
# RESULTSDIR=results \
# sbatch 1a_tblastn_and_extract.sbatch

#SBATCH -p short
#SBATCH -N 1 
#SBATCH -c 2
#SBATCH --time=12:00:00 
#SBATCH --mem=32GB                                    

# Load the necessary modules
module load ncbi-blast
module load bedtools

# Run tblastn
tblastn \
-db "$DIR/$DATABASE" \
-num_threads 2 \
-query test_query/"$QUERY" \
-out "$DATABASE"_"$QUERY".out \
-outfmt "6 qseqid sseqid bitscore qcovs evalue pident length mismatch gapopen qstart qend qseq sstart send sseq" \
-evalue 1e-5 \
-show_gis \
2> "$DATABASE"_"$QUERY".log

# Sort output
cat "$DATABASE"_"$QUERY".out \
| sort -nr -k3 \
> "$DATABASE"_"$QUERY".out.sorted

# Put the target seq hits in FASTA format
# i.e, remove gaps ("-") from the target sequence (column 15)
# then print out all of the target seqs, with the header in the form: 
# >sseqid:sstart-send	qseqid:qstart-qend	score:*,qcov:*,evalue:*,pident:*

cat "$DATABASE"_"$QUERY".out.sorted \
| awk -F" " '{gsub(/[-]/,"",$15)}1' OFS=" " \
| awk '{print ">" $2 ":" $13 "-" $14 "\t" $1 ":" $10 "-" $11 "\t" "score:" $3 ",qcov:" $4 ",eval:" $5 ",pid:" $6 "\n" $15}' \
| sed 's/ref|//g' \
| sed 's/|:/:/g' \
> "$DATABASE"_"$QUERY"_sseq.fasta

# Note: This output is the protein hit sequences, not the nucleotide sequences
# So now, extract nucleotide FASTA sequence for each hit

# Grab the headers from the amino acid FASTA file
grep '>' "$DATABASE"_"$QUERY"_sseq.fasta > "$DATABASE"_"$QUERY"_sseq_headers.txt

# Pull out all the plus strand sequences
cat "$DATABASE"_"$QUERY"_sseq_headers.txt \
| cut -f 1 \
| sed 's/>//g'  \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS=FS="\t"} {gsub(/\-/,"\t",$2)}1' \
| awk '{if ($2 < $3) print $1 "\t" $2 "\t" $3 "\t" "L1hit" "\t" "1" "\t" "+"}' \
> "$DATABASE"_"$QUERY"_plus.txt 

# Pull out all the minus strand sequences
cat "$DATABASE"_"$QUERY"_sseq_headers.txt \
| cut -f 1 \
| sed 's/>//g'  \
| sed 's/:/\t/g' \
| awk 'BEGIN {OFS=FS="\t"} {gsub(/\-/,"\t",$2)}1' \
| awk '{if ($2 > $3) print $1 "\t" $3 "\t" $2 "\t" "L1hit" "\t" "1" "\t" "-"}' \
> "$DATABASE"_"$QUERY"_minus.txt 

# Combine minus and plus strand sequences
cat "$DATABASE"_"$QUERY"_plus.txt "$DATABASE"_"$QUERY"_minus.txt  > "$DATABASE"_"$QUERY"_all.txt

# Sort and merge 
bedtools sort -i "$DATABASE"_"$QUERY"_all.txt \
| bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct > "$DATABASE"_"$QUERY"_merged.bed

# Extract FASTA from corrected merged BED file
bedtools getfasta -s -fi "$DIR/$DATABASE" -bed "$DATABASE"_"$QUERY"_merged.bed -fo "$DATABASE"_"$QUERY"_results.fasta

# Sort sequences by length
usearch -sortbylength "$DATABASE"_"$QUERY"_results.fasta -fastaout "$DATABASE"_"$QUERY"_nucl_seqs.fasta

# Move output bed and fasta to results folder
mv "$DATABASE"_"$QUERY"_nucl_seqs.fasta $RESULTSDIR
mv "$DATABASE"_"$QUERY"_merged.bed $RESULTSDIR

# Remove or move remaining files to a temporary folder
rm "$DATABASE"_"$QUERY"* 
