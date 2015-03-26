#!/usr/bin/env python
''' make_attributes.py -- converts per-residue scores to chimera attributes

    Input:
    - Tab delimited
    - Indexed by the first column
    - First row is column headers

'''

import sys
import pandas as pd
import coevo.pdb_aux.attributes as attributes
import coevo.tab_aux as tab_aux

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print >>sys.stderr, 'usage: %s scores.tab chain_id > attributes.txt' % sys.argv[0]
        sys.exit(1)
    
    res_scores_fn = sys.argv[1]
    chain_id = sys.argv[2]

    res_scores_df = tab_aux.load_flattab(res_scores_fn)

    print attributes.make_chimera_attributes(res_scores_df, chain_id)
