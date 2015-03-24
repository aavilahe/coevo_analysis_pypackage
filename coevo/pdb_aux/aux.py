#!/usr/bin/env python
''' aux.py -- auxiliary functions that help do things with Bio.PDB 

'''

import gzip
from os.path import basename
from Bio import Seq, SeqRecord, SeqUtils, PDB

from coevo.pdb_aux.distances import get_nonhet_residues

__author__ = "Aram Avila-Herrera"

def Chain_to_SeqRecord(chain):
    ''' Generates a SeqRecord from a Chain entity.

        chain: a Bio.PDB.Chain object

        Keeps only residues with blank flags (eg. no HET residues).

        Returns seqr: a Bio.SeqRecord object with a list of resnums saved in
        its letter_annotations['resnum'].

    '''

    aas = ''
    resns = list()
    for res in get_nonhet_residues(chain):
        aas += SeqUtils.seq1(res.get_resname())  # get 1-letter resname
        resns += [res.id[1]]

    seqr = SeqRecord.SeqRecord(Seq.Seq(aas), id = chain.id,
                               letter_annotations = {"resnum": resns})

    return seqr

def open_pdb(pdb_fn):
    ''' Loads structure in given pdb filename

        pdb_fn: filename as a str

        If file name ends in '.gz', attempt to open with gzip.open().

        Returns structure: a Bio.PDB.Structure object.

    '''

    if pdb_fn.endswith('.gz'):
            pdb_fh = gzip.open(pdb_fn)
    else:
            pdb_fh = open(pdb_fn)

    pdb_id = basename(pdb_fn).split('.')[0]  # set id from filename
    structure = PDB.PDBParser().get_structure(pdb_id, pdb_fh)

    return structure




