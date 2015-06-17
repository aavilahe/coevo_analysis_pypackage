#!/usr/bin/env python
''' scores.py -- load and clean raw coevo scores from many programs

    module: coevo

'''

import pandas as pd

__author__ = 'Aram.Avila-Herrera'

class Format:
    ''' Holds options for loading scores in different tabular formats:

        mfDCA
        plmDCA
        hpDCA  *formatted by hpDCA_support wrapper
        PSICOV
        CoMap
        infCalc
        distance
        CTMP
        spider (SpiderMonkey from HyPhy package)
        tab    *formatted by coevo.scores.write_tab()

    '''

    def __init__(self, prog = None, suff = '',
            offset = 0, header = None, delim = None, keep_cols = None, stat_names = None):

        # Defaults
        self.prog = prog            # tabular format name
        self.offset = offset        # alignment positions start at 0
        self.delim = delim          # None => columns are whitespace delimited
        self.keep_cols = keep_cols  # keep all the columns
        self.header = header        # no header
        self.stat_names = []        # no names
        self.preproc = None         # no preprocessing

        if prog == 'mfDCA':
            self._mfDCA_ini(suff)
        elif prog == 'plmDCA':
            self._plmDCA_ini(suff)
        elif prog == 'hpDCA':
            self._hpDCA_ini(suff)
        elif prog == 'PSICOV':
            self._PSICOV_ini(suff)
        elif prog == 'infCalc':
            self._infCalc_ini(suff)
        elif prog == 'CTMP':
            self._CTMP_ini(suff)
        elif prog == 'CoMap':
            self._CoMap_ini(suff)
        elif prog == 'spider':
            self._spider_ini(suff)
        elif prog == 'dist':
            self._distance_ini(suff)
        elif prog == 'tab':
            self._tab_ini()
    
    def _mfDCA_ini(self, suff = ''):
        ''' http://dca.ucsd.edu/DCA/DCA.html 

        '''

        self.offset = 1  # alignment positions start at 1
        self.delim  = ' '  # output columns are space-delimited
        self.header = None
        self.keep_cols = (0, 1, 2, 3)  # columns to keep
        self.stat_names = ('MIw', 'DI')
        self.stat_names = [ name + suff for name in self.stat_names ]
    
    def _plmDCA_ini(self, suff = ''):
        ''' http://plmdca.csc.kth.se

        '''

        self.offset = 1
        self.delim = ','
        self.header = None
        self.keep_cols = (0, 1, 2)
        self.stat_names = ('DIplm',)
        self.stat_names = [ name + suff for name in self.stat_names ]
    
    def _hpDCA_ini(self, suff = ''):
        ''' http://www.ploscompbiol.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pcbi.1003176.s002
        
            https://github.com/aavilahe/coevo_tools/tree/dev/wrappers

        '''
        
        self.offset = 1
        self.delim = '\t'
        self.header = None
        self.keep_cols = (0, 1, 2)
        self.stat_names = ('DI',)
        self.stat_names = [ name + suff for name in self.stat_names ]
    
    def _PSICOV_ini(self, suff = ''):
        ''' http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/

        '''

        self.offset = 1
        self.delim = ' '
        self.header = None
        self.keep_cols = (0, 1, 4)
        self.stat_names = ('PSICOV',)
        self.stat_names = [ name + suff for name in self.stat_names ]

    def _infCalc_ini(self, suff = ''):
        ''' https://github.com/aavilahe/infCalc

        '''

        self.offset = 0
        self.delim = '\t'
        self.header = 0  # first line as header
        self.keep_cols = (0, 1, 2, 3, 4, 5, 6, 7, 8)
        self.stat_names = ('Left_Entropy', 'Right_Entropy', 'Joint_Entropy', 'MI', 'VI', 'MIminh', 'MIj')
        self.stat_names = [ name + suff for name in self.stat_names ]

    def _CTMP_ini(self, suff = ''):
        ''' http://www.stat.sinica.edu.tw/chyeang/coevolution_download.zip

        '''

        self.offset = 1
        self.delim = '\t'
        self.header = 0  # first line as header
        self.keep_cols = (0, 1, 2)
        self.stat_names = ('CTMP',)
        self.stat_names = [ name + suff for name in self.stat_names ]

    def _CoMap_ini(self, suff = ''):
        ''' http://home.gna.org/comap/

        '''

        self.offset = 1
        self.delim = '\t'
        self.header = 0
        self.keep_cols = (0, 1, 5)
        self.stat_names = ('CM', 'CMP')
        self.stat_names = [ name + suff for name in self.stat_names ]
        # extract alignment positions from first column
        self.preproc = lambda df: pd.concat([df.ix[:, 0].str.extract('\[(\d+);(\d+)\]'), df.ix[:, 1:]], axis = 1) 

    def _spider_ini(self, suff = ''):
        ''' http://www.datamonkey.org/help/spidermonkey.php

            https://github.com/veg/hyphy
            https://github.com/veg/hyphy/issues/303

        '''

        self.offset = 0
        self.delim = ','
        self.header = None
        self.keep_cols = (0, 1, 2)
        self.stat_names = ('Spider',)
        self.stat_names = [ name + suff for name in self.stat_names ]

    def _distance_ini(self, suff = ''):
        self.offset = 0
        self.delim = '\t'
        self.header = None
        self.keep_cols = (0, 1, 4)
        self.stat_names = ('Dist',)
        self.stat_names = [ name + suff for name in self.stat_names ]

    def _tab_ini(self):
        self.offset = 0
        self.delim = '\t'
        self.header = 0
        self.keep_cols = None

    def load(self, fn):
        ''' load the output scores of a pairwise analysis

        '''

        df = pd.read_table(fn, sep = self.delim, header = self.header, usecols = self.keep_cols)
        if hasattr(self.preproc, '__call__'):
            df = self.preproc(df)

        idx_colnames = ['Left_Column', 'Right_Column']
        new_colnames =  idx_colnames + list(self.stat_names)
        old_colnames = df.columns[:len(new_colnames)]
        df.rename(columns = dict(zip(old_colnames, new_colnames)), inplace = True)

        df.loc[:, idx_colnames] = df.loc[:, idx_colnames].astype(int)
        df.loc[:, idx_colnames] -= self.offset  # renumber alignment columns to start at 0
        df.set_index(idx_colnames, inplace = True)

        return df

def drop_intraprotein(df, left_length):
    ''' remove intraprotein scores and renumber Right_Column to start at 0

        *assumes column numbering starts at 0*
    
    '''
    
    # pandas indices are stupid
    idx_colnames = ['Left_Column', 'Right_Column']
    interdf = df.reset_index(idx_colnames)

    # ensure Left_Column < Right_Column 
    left_vals = interdf.loc[:, idx_colnames].apply(min, axis = 1)
    right_vals = interdf.loc[:, idx_colnames].apply(max, axis = 1)
    interdf.loc[:, 'Left_Column'] = left_vals
    interdf.loc[:, 'Right_Column'] = right_vals

    # select only the interprotein rows
    interprotein_rows = ((left_vals < left_length) & (right_vals >= left_length))
    interdf = interdf[interprotein_rows].copy()
    interdf.Right_Column -= left_length  # Right_Column should count from 0

    # set the indices again
    interdf.set_index(idx_colnames, inplace = True)

    return interdf

def merge_tabs(list_of_dfs, **pd_merge_kwargs):
    ''' merge list of dataframes

    '''
    
    left_df = list_of_dfs[0]
    for right_df in list_of_dfs[1:]:
        left_df = left_df.merge(right_df, **pd_merge_kwargs)
    
    return left_df

def write_tab(df, output_fn):
    ''' save as tab delimited

    '''

    df.to_csv(output_fn, sep = '\t', float_format = '%.6f')



