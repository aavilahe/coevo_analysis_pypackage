#!/usr/bin/env python
''' fileformats.py -- format definitions for coevolution formats

    module: coevo

'''

import pandas as pd

__author__ = 'Aram.Avila-Herrera'

class Format:
    def __init__(self, prog = None, suff = '',
            offset = 0, delim = None, keep_cols = (0, 1, 2), stat_names = None):

        # defaults
        self.offset = offset  # alignment positions start at 0
        self.delim = delim  # None => whitespace delimited
        self.keep_cols = keep_cols  # keep the first three columns
        self.header = None  # no header
        if stat_names == None:
            self.stat_names = tuple([ str(sn) + suff for sn in xrange(len(keep_cols) - 2) ])
        else:
            self.stat_names = tuple([ str(sn) + suff for sn in stat_names ])
        self.preproc = None

        if prog == 'mfDCA':
            self.__mfDCA_opts__()
        elif prog == 'plmDCA':
            self.__plmDCA_opts__()
        elif prog == 'hpDCA':
            self.__hpDCA_opts__(self, suff)
        elif prog == 'PSICOV':
            self.__PSICOV_opts__()
        elif prog == 'infCalc':
            self.__infCalc_opts__()
        elif prog == 'CoMap':
            self.__CoMap_opts__()
        elif prog == 'dist':
            self.__distance_opts__()
    
    def __mfDCA_opts__(self):
        ''' http://dca.ucsd.edu/DCA/DCA.html 

        '''

        self.offset = 1  # alignment positions start at 1
        self.delim  = ' '  # output columns are space-delimited
        self.header = None
        self.keep_cols = (0, 1, 2, 3)  # columns to keep
        self.stat_names =  ('MIw', 'DI')
    
    def __plmDCA_opts__(self):
        ''' http://plmdca.csc.kth.se

        '''

        self.offset = 1
        self.delim = ','
        self.header = None
        self.keep_cols = (0, 1, 2)
        self.stat_names = ('DIplm',)
    
    def __hpDCA_opts__(self, suff = ''):
        ''' http://www.ploscompbiol.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pcbi.1003176.s002
        
            https://github.com/aavilahe/coevo_tools/tree/dev/wrappers

        '''
        
        self.offset = 1
        self.delim = '\t'
        self.header = None
        self.keep_cols = (0, 1, 2)
        self.stat_names = ('DI' + suff,)
    
    def __PSICOV_opts__(self, suff = ''):
        ''' http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/

        '''

        self.offset = 1
        self.delim = ' '
        self.header = None
        self.keep_cols = (0, 1, 4)
        self.stat_names = ('PSICOV',)

    def __infCalc_opts__(self):
        ''' https://github.com/aavilahe/infCalc

        '''

        self.offset = 0
        self.delim = '\t'
        self.header = 0  # first line as header
        self.keep_cols = (0, 1, 2, 3, 4, 5, 6, 7, 8)
        self.stat_names = ('Left_Entropy', 'Right_Entropy', 'Joint_Entropy', 'MI', 'VI', 'MIminh', 'MIj')

    def __CoMap_opts__(self, suff = ''):
        ''' http://home.gna.org/comap/

        '''

        self.offset = 1
        self.delim = '\t'
        self.header = 0
        self.keep_cols = (0, 1, 5)
        self.stat_names = ('CoMap' + suff,)
        # extract alignment positions from first column
        self.preproc = lambda df: pd.concat([df.ix[:, 0].str.extract('\[(\d+);(\d+)\]'), df.ix[:, 1:]]) 

    def __distance_opts__(self, suff = ''):
        self.offset = 0
        self.delim = '\t'
        self.header = None
        self.keep_cols = (0, 1, 4)
        self.stat_names = ('Dist' + suff,)

def load_pairwise(fn, fmt):
    ''' load the output scores of a pairwise analysis

    '''

    df = pd.read_table(fn, sep = fmt.delim, header = fmt.header, usecols = fmt.keep_cols)
    if type(fmt.preproc) == 'function':
        df = fmt.preproc(df)
    aln_column_labels = ('Left_Column', 'Right_Column')
    df.columns = aln_column_labels + fmt.stat_names
    df.loc[:, aln_column_labels] -= fmt.offset  # renumber alignment columns to start at 0
    return df

def drop_intraprotein(df, last_left_1):
    ''' remove intraprotein scores
    
    '''

    last_left_0 = last_left_1 - 1

    interprotein_rows = (df.Left_Column <= last_left_0) & (df.Right_Column > last_left_0)
    return df[interprotein_rows]

def save_tab(df, output_fn):
    df.to_csv(output_fn, sep = '\t', index = False, float_format = '%.6f')



