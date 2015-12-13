#!/usr/bin/env python
''' fasta_to_phy.py --- Converts fasta to strict phylip

    - Truncates sequence identifiers to 8 characters
    - Leaves two spaces between identifier and sequence
    - Each sequence is on its own line
    - Acts as a filter (stdin, stdout)

'''

import sys
from Bio import AlignIO

if __name__ == "__main__":
    if(len(sys.argv) > 1):
        print >>sys.stderr, "usage: %s < fasta > phy" % sys.argv[0]
        sys.exit(1)

    aln = AlignIO.read(sys.stdin, format = "fasta")
    for seq in aln:
        seq.id = seq.id[:8]
        #seq.id = seq.id[:8].ljust(8) + '  '  # explicitly truncate, pad, add two spaces
    print aln.format('phylip-sequential'),
