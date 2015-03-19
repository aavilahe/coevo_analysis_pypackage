#!/usr/bin/env python
''' concatenate_fastas.py -- concatenates two fastas horizontally

    - Sequence identifiers must be exact and unique within each fasta
    - Writes to stdout

'''

import sys
from Bio import AlignIO, SeqIO

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit('usage: %s left.fa right.fa > left_right.fa' % sys.argv[0])
        
    left_aln = AlignIO.read(sys.argv[1], format = "fasta")
    right_aln = SeqIO.to_dict(AlignIO.read(sys.argv[2], format = "fasta"))

    for left_seq in left_aln:
        if left_seq.id in right_aln:
            print (left_seq + right_aln[left_seq.id]).format('fasta'),
        else:
            print >>sys.stderr, "skipping: %s" % left_seq.id

