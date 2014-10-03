#!/usr/bin/env python
''' test coevo module

'''

import sys
import coevo  # nothing here yet
import coevo.inputformats as cf


fn1 = 'test.mfDCA'
fn2 = 'test.infCalc'
fn3 = 'test.PSICOV'

def test_load(fn, prog):
    df = cf.load_pairwise(fn, cf.Format(prog=prog))
    print df.shape
    print df[:10]
    return df

ins = list()

for fn in  [ fn1, fn2, fn3 ]:
    ins.append(test_load(fn, fn.split('.')[1]))



