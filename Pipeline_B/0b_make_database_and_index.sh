#!/bin/bash

# Index all genomes and make BLAST database
# bash 0b_make_database_and_index.sh

# Load modules if you are working on a SLURM enivronment, e.g. 
module load ncbi-blast
module load samtools

# Make each genome a BLAST database

for i in test_genome/*.fa; do makeblastdb -in $i -parse_seqids -dbtype nucl; done

# Create a fai index for each genome

for i in test_genome/*.fa; do samtools faidx $i; done
