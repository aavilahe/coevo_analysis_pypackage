#!/usr/bin/env python
''' min_dists.py -- Get min distances among chains

    When mapping distances to alignment column pairs, choose
    the closest chain.

    Eg. chain A makes contacts with a pair of identical chains C and D

'''

import sys
import pandas as pd

def load_dists(fn):
    ''' Returns pandas.DataFrame indexed by first and second column

        Expects distance name in third column
    
    '''

    df = pd.read_table(fn, header = 0, index_col = [0, 1])

    return df

def get_min(df1, df2):
    ''' Gets minimum distances for each pair of alignment columns

        pandas.concat() should work with multiple distance columns
                        but doesn't care about left/right resnums

    '''

    ## merge only works if there is a single distance column
    #dist_name = df1.columns[0]
    #dfmerge = df1.merge(df2, how = 'outer', left_index = True, right_index = True)
    #dmin = dfmerge.min(axis = 1).to_frame(dist_name)

    index_names = list(df1.index.names)
    dfcat = pd.concat([df1, df2], key = ['One', 'Two'], names = ['ChainPair'])
    dfmin = dfcat.groupby(level = index_names).min().dropna()

    return dfmin

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print >>sys.stderr, 'usage: %s dists1 dists2 > mindists'
        sys.exit(1)
    
    df1 = load_dists(sys.argv[1))
    df2 = load_dists(sys.argv[2])

    dfmin = get_min(df1, df2)
    dfmin.to_csv(sys.stdout, header = True, index = True, sep = '\t')

