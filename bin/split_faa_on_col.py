#!/usr/bin/env python
''' split_faa_on_col.py -- Splits an aligned fasta on Nth column

        col: 1,2,3,...,N,N+1,...,K
              <--Left--| |--Right-->

'''

import sys
from Bio import AlignIO

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print >>sys.stderr, 'usage: %s aln.fa nth_col left.fa right.fa' % sys.argv[0]
        sys.exit(1)
    
    aln = AlignIO.read(sys.argv[1], format = "fasta")
    nth_col = int(sys.argv[2])
    left_fn = sys.argv[3]
    right_fn = sys.argv[4]

    left_aln = aln[:,:nth_col]
    right_aln = aln[:,nth_col:]
        
    print >>open(left_fn, 'w'), left_aln.format('fasta'),
    print >>open(right_fn, 'w'), right_aln.format('fasta'),
 
