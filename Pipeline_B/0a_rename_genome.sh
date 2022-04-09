#!/bin/bash

# Rename genome file to include species name
# Example usage: GENOME=<source_genome> SPECIES=<species_name> bash 0a_rename_genome.sh

mv $GENOME "$SPECIES"_"$GENOME"
