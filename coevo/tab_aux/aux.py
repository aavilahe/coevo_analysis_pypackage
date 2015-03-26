#!/usr/bin/env python
''' aux.py -- extra utilities for dealing with tab delimited files

'''

import pandas as pd

def load_pairtab(fn):
    ''' Loads tab delimited file indexed by first two columns
        
        Returns a pandas.DataFrame

    '''

    df = pd.read_table(fn, header = 0, index_col = [0, 1])

    return df

def load_flattab(fn):
    ''' Loads tab delimited file indexed by first column

        Returns a pandas.DataFrame

    '''

    df = pd.read_table(fn, header = 0, index_col = 0)

    return df

def get_min_dists(df1, df2):
    ''' Gets minimum distances for each pair of alignment columns

        df1, df2: pandas.DataFrame indexed by alignment column numbers.
                  All DataFrame columns are expected to be distances

        pandas.concat() should work with multiple distance columns
                        but doesn't care about index names

    '''

    ## merge only works if there is a single distance column
    #dist_name = df1.columns[0]
    #dfmerge = df1.merge(df2, how = 'outer', left_index = True, right_index = True)
    #dmin = dfmerge.min(axis = 1).to_frame(dist_name)

    index_names = list(df1.index.names)
    dfcat = pd.concat([df1, df2], keys = ['One', 'Two'], names = ['ChainPair'])
    dfmin = dfcat.groupby(level = index_names).min().dropna()

    return dfmin

def convert_col(target_df, map_df, from_col, to_col):
    ''' Converts from_col to to_col in target_df using map_df

        target_df: pandas.DataFrame containing from_col in column
        map_df: pandas.DataFrame mapping from_col to to_col
        from_col, to_col: column names to convert

        Returns target_df with from_col replaced with to_col

    '''

    target_dfc = target_df.merge(map_df, how = 'inner', on = from_col
                                 ).drop([from_col], axis = 1)

    return target_dfc

