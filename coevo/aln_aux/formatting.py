#!/usr/bin/env python
''' formatting.py -- utilities for formatting sequence alignments

'''

def make_strict_phylip_id_map(seqids):
    ''' Maps strict phylip ids to original ids
    
        seqids: list of seqid strings

        Strict phylip ids are 10 chars long. The last two characters
        are always spaces.

        eg. "mosquito_virus_[strain:36f]" --> "32_mosqu  "

        Returns a generator of tuples: ((old_id, new_id), ...)

    '''

    return (((seqid, str(i) + '_' +  seqid)[:8].ljust(8) + '  ')
            for i, seqid
            in enumerate(seqids)
            )

def replace_ids(aln, id_map):
    ''' Replaces ids in aln using id_map

        aln: Bio.Align.MultipleSeqAlignment object
        id_map: dict from old to new ids

    '''

    for seqr in aln:
        seqr.id = id_map[seqr.id]
    
    return aln

