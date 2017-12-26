#!/bin/bash

# $1=input fasta, $2=window size, $3=step size, $4=k size, $5=label

fragments="$(./fragmentize $1 $2 $3)"

./fasta_to_features.py "$1_frag.fa" $fragments $4 "$1_raw.csv"
./feat_prep "$1_raw.csv" "$1_prep.csv" $5

