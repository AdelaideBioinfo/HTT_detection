#!/bin/bash

# Append species name to all sequence headers
#
# Example usage:
#
# SPECIES=YarrowiaLipolytica \
# ELEMENT=L1 \
# RESULTSDIR=results \
# bash 1d_append_name_to_headers.sh

cat $RESULTSDIR/"$SPECIES"_"$ELEMENT"_combined.fasta \
| sed 's/^>/>'$SPECIES'_/g' \
| awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' \
> $RESULTSDIR/"$SPECIES"_"$ELEMENT"_final.fasta

