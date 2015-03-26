#!/usr/bin/python env
''' aux.py -- miscellaneous utility functions

'''

import sys

import tempfile
from Bio import SeqIO, Align
from copy import deepcopy

def pop_row(aln, seqid):
    ''' Pop a row from an alignment by sequence id

        aln: a Bio.Align.MultipleSeqAlignment object
        seqid: id of Bio.SeqRecord.SeqRecord to pop from aln

        Returns a tuple containing the popped SeqRecord and a
        copy of aln without seqid's SeqRecord.

    '''

    aln_d = SeqIO.to_dict(aln)
    seq = aln_d[seqid]
    del aln_d[seqid]
    aln = Align.MultipleSeqAlignment(aln_d.itervalues())

    return seq, aln

def annotate_positions(seqrec):
    ''' Adds position numbering to a SeqRecord
        
    '''

    pos_list = list()
    pos = 0
    for aa in seqrec.seq:
        if aa == '-':
            pos_list += ['-']
        else:
            pos_list += [pos]
            pos += 1

    seqrec.letter_annotations['pos'] = pos_list

    return seqrec

def ungap_SeqRecord(seqr):
    ''' Returns an ungapped copy of a SeqRecord

    '''

    ug_seqr = deepcopy(seqr)
    ug_seqr.seq = ug_seqr.seq.ungap('-')

    return ug_seqr

def make_tmp_fa(seqrec):
    ''' Writes SeqRecord to temporary file in fasta format

        Returns handle to temp file

    '''

    tmp_fh = tempfile.NamedTemporaryFile(suffix='.fa', delete=False)
    SeqIO.write([seqrec], tmp_fh, 'fasta')
    print >>sys.stderr, 'Created temporary fasta "%s"' % tmp_fh.name
    tmp_fh.close()

    return tmp_fh

