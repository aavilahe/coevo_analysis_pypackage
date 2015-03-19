#!/usr/bin/env python
''' concatenate_fastas.py -- concatenates two fastas horizontally

    - Sequence identifiers must be exact and unique within each fasta
    - Writes to stdout

'''

import sys
from Bio import AlignIO

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit('usage: %s left.fa right.fa > left_right.fa' % sys.argv[0])
        
    left_aln = AlignIO.read(sys.argv[1], format = "fasta")
    right_aln = AlignIO.read(sys.argv[2], format = "fasta")

    right_map = dict(((seq.id, seqn) for seqn, seq in enumerate(right_aln)))

    for left_seq in left_aln:
        if left_seq.id in right_map:
            print (left_seq + right_aln[right_map[left_seq.id]]).format('fasta'),
        else:
            print >>sys.stderr, "skipping: %s" % left_seq.id

