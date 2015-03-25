#!/usr/bin/env python
''' distances.py -- for working with distances in pdb files

'''

from numpy.linalg import norm
from itertools import product

__author__ = 'Aram Avila-Herrera'

def get_CB_coord(residue):
    ''' Get beta carbon coordinates from residue

        residue: a Bio.PDB.Residue
        
        If residue is a glycine, returns alpha carbon coordinates
        Returns a list containing one coordinate as a numpy.ndarray of float32

    '''

    if residue.get_resname() == 'GLY':
        return [residue['CA'].coord]
    else:
        return [residue['CB'].coord]

def get_nonH_coords(residue):
        ''' Get coordinates from non-hydrogen atoms in residue
            residue: a Bio.PDB.Residue

            Returns a list of coordinates, each a numpy.ndarray

        '''

        return [atom.coord
                for atom
                in residue.get_list()
                if not atom.get_id().startswith('H')
                ]

def get_allatom_coords(residue):
    ''' Get coordinates from all atoms in residue
        
        residue: a Bio.PDB.Residue

        Returns a list of coordinates, each a numpy.ndarray

    '''

    return [atom.coord
            for atom
            in residue.get_list()
            ]

def has_structure_carbon(residue):
    ''' Returns True if residue has a CB (CA if residue is a glycine)
        
    '''

    if 'CB' in residue:
        return True
    if residue.get_resname() == 'GLY' and 'CA' in residue:
        return True

    return False

def get_nonhet_residues(chain):
    ''' Get residues in a Bio.PDB.Chain entity with blank HET flag

        Returns a generator over those residues

    '''

    residues = chain.get_residues()
    nonhet_res = (res for res in residues if res.id[0] == ' ')

    return nonhet_res

def calc_residue_distance(res_a, res_b, get_coords):
    ''' Get euclidean distance between residues

        res_a, res_b: Bio.PDB.Residue entities
        get_coords: a function that returns atom coordinates for a residue
                    get_coords() should return a list of numpy.ndarrays

        Returns a numpy.float32
        
    '''
    
    coords_a = get_coords(res_a)
    coords_b = get_coords(res_b)

    atom_distances = (
                      norm(X1 - X2)
                      for (X1, X2)
                      in product(coords_a, coords_b)
                      )
    
    # We want the smallest atom-atom distance between residues
    return min(atom_distances)

def choose_get_coords(dist_atoms):
    ''' Choose get_coords function from dist_atoms

        dist_atoms: a string 'Cb', 'NoH', or 'Any'

        Returns a get_coords function. Defaults to get_CB_coord()
        
        TODO: make these part of a class or something...

    '''
    
    if dist_atoms == 'NoH':
        return get_nonH_coords
    if dist_atoms == 'Any':
        return get_allatom_coords
    return get_CB_coord

