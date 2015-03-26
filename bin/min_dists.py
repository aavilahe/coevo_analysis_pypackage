#!/usr/bin/env python
''' min_dists.py -- Get min distances among chains

    When mapping distances to alignment column pairs, choose
    the closest chain.

    Eg. chain A makes contacts with a pair of identical chains C and D

'''

import sys
import pandas as pd
import coevo.pdb_aux.distances as distances

def load_dists(fn):
    ''' Returns pandas.DataFrame indexed by first and second column
        
        Expects distances are in third column. Drops other columns

    '''
    
    return distances.load_pairtab(fn).ix[:, :3]

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >>sys.stderr, 'usage: %s dists1 dists2 > mindists' % sys.argv[0]
        sys.exit(1)
    
    df1 = load_dists(sys.argv[1])
    df2 = load_dists(sys.argv[2])
    
    dfmin = distances.get_min_dists(df1, df2)
    dfmin.to_csv(sys.stdout, header = True, index = True,
                 sep = '\t', float_format = '%.6f'
                 )

