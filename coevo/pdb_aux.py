#!/usr/bin/env python
''' pdb_aux.py -- auxiliary functions that help do things with Bio.PDB 

'''

from Bio import Seq, SeqRecord, SeqUtils

__author__ = "Aram Avila-Herrera"

def Chain_to_SeqRecord(chain):
    ''' Generates a SeqRecord from a Chain entity.

        chain: a Bio.PDB.Chain object

        Keeps only residues with blank flags (eg. no HET residues).

        Returns seqr: a Bio.SeqRecord object with a list of resnums saved in
        its letter_annotations['resnum'].

    '''

    # aa_l, resns = zip(*[
    #                       (SeqUtils.seq1(res.resname), res.id[1])
    #                       for res
    #                       in chain.get_residues()
    #                       if res.id[0] == ' '
    #                     ])
    # aas = ''.join(aa_l)

    residues = chain.get_residues()

    aas = ''
    resns = list()
    for res in residues:
        # res.id is currently a 3-tuple: (het, resnum, icode)
        if res.id[0] == ' ':
            # no HET flag
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




