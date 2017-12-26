#!/bin/bash

# $1=input fasta, $2=window size, $3=step size, $4=k size, $5=label, $6 typeI error, $7 typeII error

rm hmmout.txt

echo Fragmentizing input...
fragments="$(./fragmentize $1 $2 $3)"

echo "Preparing $4-mer table..."
./fasta_to_features.py "$1_frag.fa" $fragments $4 kmer_raw.csv
./feat_prep kmer_raw.csv kmer_prep.csv $5

echo Predicting class labels...
/home/zhouw/R/bin/Rscript predict_label.R

echo Calling class labels using HMM...
./obs_states kmer_raw.csv predoutput.txt $6 $7
sort -t ',' -k10 -r -n HMM_states.csv > HMM_sorted.csv
/home/zhouw/R/bin/Rscript predict_HMM.R
sort hmmout.txt > hmmout_sorted.txt

echo Extracting viral sequences...
./extract_viral_fasta "$1_frag.fa"
mv viral_sequences.fa viral_sequences-$1.fa

