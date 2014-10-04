#!/usr/bin/env python
''' test coevo module

'''

import sys
import coevo  # nothing here yet
import coevo.io as cio


fn1 = 'test.mfDCA'
fn2 = 'test.infCalc'
fn3 = 'test.PSICOV'

def dfstat(df):
    print df.shape
    print df[:10]

def test_load(fn, prog):
    df = cio.load_pairwise(fn, cio.Format(prog=prog))
    dfstat(df)
    return df

print '\n\n ..... test load \n\n'
ins = list()
for fn in  [ fn1, fn2, fn3 ]:
    ins.append(test_load(fn, fn.split('.')[1]))

print '\n\n ..... drop intra \n\n'

for i in xrange(len(ins)):
    df = ins[i]
    print '...'
    if i != 1:
        df = cio.drop_intraprotein(df, 310)
    dfstat(df)

for i,df in enumerate(ins):
    cio.save_tab(df, str(i) + '.tab')

