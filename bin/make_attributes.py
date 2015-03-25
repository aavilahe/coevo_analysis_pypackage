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

def load_scores(fn):
    ''' Returns pandas.DataFrame indexed by first column

        Column names defined in first row
    
    '''

    df = pd.read_table(fn, header = 0, index_col = 0)

    return df


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print >>sys.stderr, 'usage: %s scores.tab chain_id > attributes.txt' % sys.argv[0]
        sys.exit(1)
    
    res_scores_fn = open(sys.argv[1], 'r')
    chain_id = sys.argv[2]

    res_scores_df = load_scores(res_scores_fn)

    print attributes.make_chimera_attributes(res_scores_df, chain_id)
