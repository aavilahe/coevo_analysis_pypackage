#!/usr/bin/env python
''' map_column_to_resnum.py -- map alignment column to resnum in pdb chain

    Example: chain A in 3DGE starts at 246, but alignment covers discontinuous domains

'''

import sys
import getopt

from Bio import SeqRecord, Align, AlignIO

import coevo.pdb_aux as pdb_aux
import coevo.aln_aux as aln_aux
import coevo.aln_aux.wrappers as aln_wraps

__author__ = "Aram Avila-Herrera"

PROFILE_ALIGNER_CMD = 'fmuscle -profile -in1 %s -in2 %s -out /dev/stdout'
PAIR_ALIGNER_CMD = 'needle -auto -asequence %s -bsequence %s -stdout -aformat3 fasta'

def parse_cmd_line(args):
    ''' Parse command line for options and return a dict.

        uses getopt for compatibility

    '''

    optlist, args = getopt.getopt(args, 'hr:i',
                                  ['help',
                                   'refid=',
                                   'int_aln',
                                   'user_aln'
                                   ]
                                  )
    options = dict()
    usage = (
             'usage: %s [ options ] chain_id pdb_file aln.fa\n\n'
             'Options:\n'
             '-h, --help           print this help and exit\n'
             '-r, --refid REFID    map using given reference seq id in alignment\n'
             '-i, --int_aln        use internal global aligner\n'
             '                     (ignored if -r|--refid is not specified)\n'
             '--user_aln           use user defined aligners\n'
             '                     See: PROFILE_ALIGNER_CMD and PAIR_ALIGNER_CMD\n\n'

             ) % sys.argv[0]

    # Read command line
    for opt, val in optlist:
        if opt in ('-h', '--help'):
            sys.exit(usage)
        if opt in ('-r', '--refid'):
            options['refid'] = val
        if opt in ('--user_aln'):
            options['user_aln'] = True
        if opt in ('-i', '--int_aln'):
            options['int_aln'] = True

    if len(args) != 3:
        print >>sys.stderr, 'wrong number of arguments'
        sys.exit(usage)

    options['chain_id'] = args[0]
    options['pdb_file'] = args[1]
    options['aln.fa'] = args[2]

    return options

def mangle_id(seq_id):
    ''' mangle sequence id to avoid conflicts

    '''

    return '__XX__chain[%s]__XX__' % seq_id

def col_to_resn(aln, chain_seqr, aln_ref, chain_ref_id):
    ''' map columns in aln to resnums in chain_seqr

        aln: original input alignment
        chain_seqr: chain sequence annotated with resnums
        aln_ref: aln aligned to chain_seqr (profile or pairwise to ref)
        chain_ref_id: mangled chain sequence id in aln_ref

        Returns
        col_resn_aa: list of (column number, resnum, amino acid) tuples

        Mapping is done in three steps:

        1. map aln column to aln_ref column
        2. map aln_ref column to chain_ref position
        3. map chain_ref position to resnum (chain_ref position = chain_seqr position)

    '''

    chain_id = chain_seqr.id

    # pop aligned chain sequence from reference alignment
    chain_ref, aln_ref = aln_aux.pop_row(aln_ref, chain_ref_id)
    chain_ref = aln_aux.annotate_positions(chain_ref)

    # sort aln_ref and aln rows by sequence id (to compare columns)
    aln_ref.sort()
    aln.sort()

    # 1. map aln column to aln_ref columns
    aln_len = aln.get_alignment_length()
    aln_ref_len = aln_ref.get_alignment_length()
    js = 0  # first column in aln_ref to start checking
    i_j = list()
    for i in xrange(aln_len):
        if str(set(aln[:, i])) != '-':
            # all-gap columns in original alignment result in ambiguous mapping
            # so we skip them
            for j in xrange(js, aln_ref_len):
                if aln_ref[:, j] == aln[:, i]:
                    i_j += [(i, j)]
                    js = j + 1
                    break

    if len(i_j) != aln_len:
        print >>sys.stderr, 'Warning: some original alignment columns were skipped while mapping!'

    # 2. map aln_ref column to exchain position
    # 3. map chain_seqr index to resnum
    col_resn_aa = list()

    for (i, j) in i_j:
        pos = chain_ref.letter_annotations['pos'][j]
        if pos != '-':
            resn = chain_seqr.letter_annotations['resnum'][pos]
            col_resn_aa += [(i, resn, chain_seqr[pos])]

    return col_resn_aa


if __name__ == "__main__":
    options = parse_cmd_line(sys.argv[1:])
    aln = AlignIO.read(options['aln.fa'], format = "fasta")
    structure = pdb_aux.open_pdb(options['pdb_file'])
    model = structure[0]
    chain = model[options['chain_id']]
    chain_seqr = pdb_aux.Chain_to_SeqRecord(chain)
    chain_seqr_m = SeqRecord.SeqRecord(chain_seqr.seq,
                                       id = mangle_id(chain_seqr.id),
                                       name = '', description = ''
                                       )

    if 'refid' in options:
        # align chain to reference in alignment
        aln_l = [ seqr for seqr in aln if seqr.id == options['refid'] ]
        if len(aln_l) < 1:
            print >>sys.stderr, "refid <%s> not in reference alignment %s" % \
                                (options['refid'], options['aln.fa'])
            sys.exit(1)
        ref_seqr = aln_l[0]
        if 'user_aln' in options:
            ex_aligner = aln_wraps.make_external_aligner(PAIR_ALIGNER_CMD)
        elif 'int_aln' in options:
            ex_aligner = None
        else:
            ex_aligner = aln_wraps.needle_align
        aln_ref = aln_wraps.pair_align_SeqRecords(chain_seqr_m, ref_seqr, ex_aligner)
        col_resn_aa = col_to_resn(Align.MultipleSeqAlignment([ref_seqr]),
                                  chain_seqr, aln_ref, chain_seqr_m.id
                                  )
    else:
        # align chain to alignment
        if 'user_aln' in options:
            ex_aligner = aln_wraps.make_external_aligner(PROFILE_ALIGNER_CMD)
        else:
            ex_aligner = aln_wraps.muscle_profile_align
        aln_ref = aln_wraps.profile_align_SeqRecord_to_fa(chain_seqr_m,
                                                    options['aln.fa'],
                                                    ex_aligner
                                                    )
        col_resn_aa = col_to_resn(aln, chain_seqr, aln_ref, chain_seqr_m.id)

    print '\t'.join(['Column', 'resn', 'AA'])
    for c_r_a in col_resn_aa:
        print '\t'.join(map(str, c_r_a))

