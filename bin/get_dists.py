#!/usr/bin/env python
''' get_dists.py -- Gets intra- or interchain residue-residue distances

    User can specify distance type:
    - beta carbon (default)
    - nearest non-hydrogen atom
    - nearest atom

    Outputs resnums, resnames, distances

    Can output alignment columns if mappings from resnums to alignment columns
    are specified

'''

import sys
import getopt
from itertools import product, combinations
import pandas as pd
from Bio import SeqUtils

import coevo.pdb_aux as pdb_aux
import coevo.tab_aux as tab_aux
import coevo.pdb_aux.distances as distances

__author__ = 'Aram Avila-Herrera'

def parse_cmd_line(args):
    ''' Parse command line or config file for options and return a dict.

        Uses getopt for compatibility

    '''

    optlist, args = getopt.getopt(args, 'hc:d:',
                                  [
                                   'help',
                                   'chainL=',
                                   'chainR=',
                                   'mapL=',
                                   'mapR=',
                                   'dist_atoms='
                                    ]
                                  )
    options = dict()
    usage = (
             'usage: %s [ options ] pdb_file\n\n'
             'Options:\n'
             '\t-h, --help                       show this help message and exit\n'
             '\t-c, --chainL <chain id>          (required) specify single or left chain id\n'
             '\t--chainR <chain id>              specify right chain id\n'
             '\t-d, --dist_atoms Cb | NoH | Any  specify beta carbon (default), non-hydrogen, or any atom distances\n'
             '\t--mapL <resn2col_left>           specify a mapping from resnum to alignment columns for left chain\n'
             '\t--mapR <resn2col_right>          specify a mapping from resnum to alignment columns for right chain\n'
             ) % sys.argv[0]

    # Set defaults
    options['distance'] = 'Cb'

    # Assign options
    for opt, val in optlist:
        if opt in ('-h', '--help'):
            sys.exit(usage)
        if opt in ('-c', '--chainL'):
            options['chainL'] = val
        if opt == '--chainR':
            options['chainR'] = val
        if opt in ('-d', '--dist_atoms'):
            options['dist_atoms'] = val
        if opt in ('--mapL'):
            options['mapL'] = val
        if opt in ('--mapR'):
            options['mapR'] = val

    # Check arguments and required options
    if len(args) < 1:
        err = 'Error: specify a pdb_file'
        sys.exit(err + '\n' + usage)
    options['pdb_file'] = args[0]

    if 'chainL' not in options:
        err = 'Error: specify a single or left chain'
        sys.exit(err + '\n' + usage)

    return options

def get_residues(chain):
    ''' Returns generator over structural residues in chain

        chain: Bio.PDB.Chain entity

    '''

    return (res
            for res
            in distances.get_nonhet_residues(chain)
            if distances.has_structure_carbon(res)
            )

def make_interchain_pairs(chain_a, chain_b):
    ''' Return a generator over pairs of residues in two chains.

        chain_a, chain_b: Bio.PDB.Chain entity

    '''

    res_a_g = get_residues(chain_a)
    res_b_g = get_residues(chain_b)
    return product(res_a_g, res_b_g)

def make_intrachain_pairs(chain):
    ''' Return a generator over pairs of residues in a chain

        chain: Bio.PDB.Chain entity

    '''

    return combinations(get_residues(chain), r = 2)

def get_distances(res_pairs, get_coords):
    ''' Get distances for all pairs of residues between two chains

        res_pairs: generator over tuples ((res_a, res_b), ...)
        get_coords: function to get residue coordinates

        Returns a list over 5-tuples: [(resn_a, resn_b, aa_a, aa_b, dist), ...]

    '''

    return [
            (res_a.id[1], res_b.id[1],
             distances.calc_residue_distance(res_a, res_b, get_coords),
             SeqUtils.seq1(res_a.resname), SeqUtils.seq1(res_b.resname)
             )
            for (res_a, res_b)
            in res_pairs
            ]


if __name__ == "__main__":
    options = parse_cmd_line(sys.argv[1:])
    structure = pdb_aux.open_pdb(options['pdb_file'])
    model = structure[0]
    get_coords = distances.choose_get_coords(options['dist_atoms'])

    if 'chainR' in options:
        res_pairs = make_interchain_pairs(model[options['chainL']],
                                         model[options['chainR']]
                                         )
    else:
        res_pairs = make_intrachain_pairs(model[options['chainL']])

    dists = get_distances(res_pairs, get_coords)
    dists_df = pd.DataFrame(dists,
                            columns = ['Left_resn', 'Right_resn',
                                       'Distance',
                                       'Left_AA', 'Right_AA'
                                       ]
                            )

    if 'mapL' in options:
        lmap = tab_aux.load_map(options['mapL'], 'Left')  # Assumes "Column" is 1st column in header
        dists_df = tab_aux.convert_col(dists_df, lmap,
                                       from_col = 'Left_resn',
                                       to_col = 'Left_Column'
                                       )
    if 'mapR' in options:
        rmap = tab_aux.load_map(options['mapR'], 'Right')
        dists_df = tab_aux.convert_col(dists_df, rmap,
                                       from_col = 'Right_resn',
                                       to_col = 'Right_Column'
                                       )

    indices = [col
               for col
               in dists_df.columns
               if col in ('Left_resn', 'Right_resn',
                          'Left_Column', 'Right_Column'
                          )
               ]
    dists_df.set_index(indices, inplace = True)
    dists_df.to_csv(sys.stdout, sep = '\t', header = True, float_format = '%.6f')

