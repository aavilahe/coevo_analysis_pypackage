#!/usr/bin/env python
''' test coevo module

'''

import sys
import coevo  # nothing here yet
import coevo.scores as cs


def diag(df, astr=''):
    print astr
    print df[:5]
    print df.shape
    print '\n\n'

fnlist = [
    'test.mfDCA',
    'test.infCalc',
    'test.PSICOV'
    ]

tabdf_dict = dict()
for fn in fnlist:
    print "Reading %s" % fn
    prog = fn.split('.')[-1]
    tabdf = cs.Format(prog = prog).load(fn)
    tabdf_dict[prog] = tabdf
    diag(tabdf, '%s : ' % fn)
    
for prog, df in tabdf_dict.iteritems():
    print "Dropping intraprotein in %s" % prog
    if prog != 'infCalc':
        tabdf_dict[prog] = cs.drop_intraprotein(df, left_length = 310)
    diag(tabdf_dict[prog], '%s : ' % prog)

bigdf = cs.merge_tabs(tabdf_dict.values(), left_index = True, right_index = True, how = 'outer')
diag(bigdf, 'merged:')

cs.write_tab(bigdf, 'big_test.tab')






