#!/usr/bin/env python
''' convert_resnums_to_columns.py -- Converts resnums in get_dists.py output

'''

import sys

import pandas as pd
import coevo.pdb_aux.distances as distances

def load_map(fn, left_or_right):
    ''' Loads column-to-resn map and adds prefix

        Drops extra columns

    '''

    df = distances.load_pairtab(fn).reset_index().ix[:, :2]
    df.rename(columns = {'resn': left_or_right + '_resn',
                         'Column': left_or_right + '_Column'
                         },
              inplace = True
              )
    return df

def convert_resns(dist_df, leftmap_df, rightmap_df):
    ''' Converts resnum in dist_df to alignment column numbers

    '''
    
    dist_dfm = dist_df.merge(leftmap_df,
                             how = 'inner',
                             on = 'Left_resn'
                             ).merge(rightmap_df,
                                     how = 'inner',
                                     on = 'Right_resn'
                                     )
    dist_dfm.drop(['Left_resn', 'Right_resn'], axis = 1, inplace = True)
    dist_dfm.set_index(['Left_Column', 'Right_Column'], inplace = True)

    return dist_dfm

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print >>sys.stderr, 'usage: %s resn_dists left_map right_map > coln_dists' % sys.argv[0]
        sys.exit(1)
    
    dist_df = distances.load_pairtab(sys.argv[1]).reset_index()    
    lmap_df = load_map(sys.argv[2], 'Left')
    rmap_df = load_map(sys.argv[3], 'Right')

    df = convert_resns(dist_df, lmap_df, rmap_df)
    df.to_csv(sys.stdout, header = True, index = True,
              sep = '\t', float_format = '%.6f'
              )

