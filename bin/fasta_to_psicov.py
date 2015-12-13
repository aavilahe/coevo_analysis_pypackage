#!/usr/bin/env python
''' fasta_to_psicov.py --- Converts fasta to PSICOV input format

    - Removes sequence identifiers
    - PSICOV format: one sequence per line, no whitespace
    - Acts as filter (stdin, stdout)

'''

import sys
from Bio import AlignIO

if __name__ == "__main__":
    if(len(sys.argv) > 1):
        print >>sys.stderr, "usage: %s < fasta > psicov" % sys.argv[0]
        sys.exit(1)
    aln = AlignIO.read(sys.stdin, format = "fasta")
    aln = '\n'.join([str(seq.seq) for seq in aln])
    print aln,
