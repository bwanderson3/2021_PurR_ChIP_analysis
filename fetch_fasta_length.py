#! /usr/bin/env python3

from Bio import SeqIO

seqs = [i for i in SeqIO.parse(input("Type the full path to your fasta file here: "), 'fasta')]

for seq in seqs:
    print(len(seq))