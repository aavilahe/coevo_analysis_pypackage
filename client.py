#!/usr/bin/env python
''' test coevo module

'''

import sys
import coevo  # nothing here yet
import coevo.inputformats as cf


fn1 = 'test.mfDCA'
fn2 = 'test.infCalc'

fmt1 = cf.Format(prog='mfDCA')
fmt2 = cf.Format(prog='infCalc')

in1 = cf.load_pairwise(fn1, fmt1)
in2 = cf.load_pairwise(fn2, fmt2)






