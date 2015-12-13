#!/usr/bin/env python
''' make_attributes.py --- converts per-residue scores to chimera attributes

    Converts scores for residues on a single chain to chimera attributes.

    Input is:
    - Tab delimited
    - Indexed by the first column
    - First row is column headers

'''

import sys
import getopt
import pandas as pd
import coevo.pdb_aux.attributes as attributes
import coevo.tab_aux as tab_aux

def parse_cmd_line(args):
    ''' Parse command line or config file for options and return a dict.

        Uses getopt for compatibility

    '''

    optlist, args = getopt.getopt(args, 'hc:m:',
                                  [
                                   'help',
                                   'chain_id=',
                                   'map=',
                                    ]
                                  )
    options = dict()
    usage = (
             'usage: %s [ options ] scores_tab\n\n'
             'Options:\n'
             '\t-h, --help                    show this help message and exit\n'
             '\t-c, --chain_id <chain id>     specify single or chain id\n'
             '\t--map <resn2col>              specify a mapping from resnum to alignment columns for chain\n'
             ) % sys.argv[0]

    # Assign options
    for opt, val in optlist:
        if opt in ('-h', '--help'):
            sys.exit(usage)
        if opt in ('-c', '--chain_id'):
            options['chain_id'] = val
        if opt in ('--map'):
            options['map'] = val

    # Set some defaults
    if 'chain_id' not in options:
        options['chain_id'] = ''

    # Check arguments and required options
    if len(args) < 1:
        err = 'Error: specify a tab delimited file containing scores'
        sys.exit(err + '\n' + usage)
    options['scores_tab'] = args[0]

    return options

if __name__ == "__main__":
    options = parse_cmd_line(sys.argv[1:])

    res_scores_fn = options['scores_tab']
    res_scores_df = tab_aux.load_flattab(res_scores_fn)

    if 'map' in options:
        res_scores_df.reset_index(inplace = True)
        map_df = tab_aux.load_map(options['map'])
        res_scores_df = tab_aux.convert_col(res_scores_df, map_df,
                                            from_col = 'Column',
                                            to_col = 'resn'
                                            )
        res_scores_df.set_index(['resn'], inplace = True)

    # if chain_id = '' if not user-specified
    print attributes.make_chimera_attributes(res_scores_df, options['chain_id'])
