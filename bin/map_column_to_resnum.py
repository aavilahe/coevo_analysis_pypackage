#!/usr/bin/env python
''' map_column_to_resnum.py -- map alignment column to resnum in pdb chain

    Example: chain A in 3DGE starts at 246, but alignment covers discontinuous domains

'''

import sys
import getopt
import gzip
import tempfile
from os import remove
from os.path import basename
from subprocess import *
from StringIO import StringIO

from Bio import Seq, SeqRecord, SeqIO, Align, AlignIO, PDB, SeqUtils
from Bio.Align.Applications import MuscleCommandline
from Bio.Emboss.Applications import NeedleCommandline

__author__ = "Aram Avila-Herrera"

PROFILE_ALIGNER_CMD = 'fmuscle -profile -in1 %s -in2 %s -out /dev/stdout'
PAIR_ALIGNER_CMD = 'needle -auto -asequence %s -bsequence %s -stdout -aformat3 fasta'

def parse_cmd_line(args):
    ''' Parse command line for options and return a dict.

        uses getopt for compatibility

    '''

    optlist, args = getopt.getopt(args, 'hr:', ['help', 'refid=', 'ex_aln'])
    options = dict()
    usage = (
             'usage: %s [ options ] chain_id pdb_file aln.fa\n\n'
             'Options:\n'
             '-h, --help           print this help and exit\n'
             '-r, --refid REFID    map using given reference seq id in alignment\n'
             '--ex_aln             use external aligners\n\n'
             ) % sys.argv[0]

    # Read command line
    for opt, val in optlist:
        if opt in ('-h', '--help'):
            sys.exit(usage)
        if opt in ('-r', '--refid'):
            options['refid'] = val
        if opt in ('--ex_aln'):
            options['ex_aln'] = True

    if len(args) != 3:
        print >>sys.stderr, 'wrong number of arguments'
        sys.exit(usage)

    options['chain_id'] = args[0]
    options['pdb_file'] = args[1]
    options['aln.fa'] = args[2]

    return options

def open_model(pdb_fn):
    ''' opens first model in given pdb file name

        if file name ends in .gz, will attempt to open as a gzipped file

    '''

    if pdb_fn.endswith('.gz'):
            pdb_fh = gzip.open(pdb_fn)
    else:
            pdb_fh = open(pdb_fn)

    pdb_id = basename(pdb_fn).split('.')[0]
    structure = PDB.PDBParser().get_structure(pdb_id, pdb_fh)
    model = structure[0] # assume 1 model in pdb file

    return model

def chain_to_SeqRecord(chain, chain_id):
    ''' generate a SeqRecord from chain object and chain id

        keeps only residues with no flags (eg. no HET residues)
        resnum saved in letter_annotations['resnum']

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

    seqrec = SeqRecord.SeqRecord(Seq.Seq(aas), id = chain_id,
                           letter_annotations = {"resnum": resns})

    return seqrec

# aln
def make_ex_aligner(EX_ALN_CMD):
    ''' returns a function to align to fasta files with an external aligner

        EX_ALN_CMD: command line containing two '%s' formatting operators, eg.
                    'muscle -profile -in1 %s -in2 %s -out /dev/stdout' or
                    'needle -auto -asequence %s -bsequence %s -stdout -aformat3 fasta'

    '''

    def ex_aligner(fa1, fa2):
        ''' aligns sequences in two fasta files with an external aligner

            fa1, fa2: filenames of fasta files to align

            returns a MultipleSeqAlignment object

        '''

        cmd = (EX_ALN_CMD % (fa1, fa2)).split()
        cmd_sout = Popen(cmd, stdout=PIPE).communicate()[0]  # watch out for large alignments
        exaln = AlignIO.read(StringIO(cmd_sout), format = "fasta")

        return exaln

    return ex_aligner

def muscle_profile_align(fa1, fa2):
    ''' uses muscle to align two fastas

        fa1 and fa2 must exist on disk when command is called

        returns a MultipleSeqAlignment object

    '''

    muscle_cmd = MuscleCommandline(input = fa1, in2 = fa2,
                                   profile = True
                                   )
    exaln = AlignIO.read(StringIO(muscle_cmd()[0]), format = "fasta")

    return exaln

def needle_align(fa1, fa2):
    ''' uses needle to align two fastas

        fa1 and fa2 must exist on disk when command is called

        default penalties:
            gapopen = 10.0
            gapextend = 0.5

        returns a MultipleSeqAlignment object

    '''

    needle_cmd = NeedleCommandline(asequence = fa1, bsequence = fa2,
                                   outfile='/dev/stdout', aformat = 'fasta',
                                   gapopen = 10.0, gapextend = 0.5
                                   )
    exaln = AlignIO.read(StringIO(needle_cmd()[0]), format = 'fasta')

    return exaln

def profile_align(seqr, aln_fn, ex_aligner = muscle_profile_align):
    ''' align sequence to alignment

    '''

    tmp_fa = make_tmp_fa(seqr.format('fasta'))
    exaln = ex_aligner(tmp_fa.name, aln_fn)
    remove(tmp_fa.name)

    return exaln

def pair_align(seqr, ref_seqr, ex_aligner = needle_align):
    ''' align sequence to reference sequence in alignment

    '''

    tmp_fa = make_tmp_fa(seqr.format('fasta'))
    tmp_ref_fa = make_tmp_fa(SeqRecord.SeqRecord(ref_seqr.seq.ungap('-'),
                             id = ref_seqr.id, description = '').format('fasta'))
    exaln = ex_aligner(tmp_fa.name, tmp_ref_fa.name)
    remove(tmp_fa.name)
    remove(tmp_ref_fa.name)

    return exaln

def mangle_id(seq_id):
    ''' mangle sequence id to avoid conflicts

    '''

    return '__XX__chain[%s]__XX__' % seq_id

def make_tmp_fa(fasta_str):
    ''' writes fasta formatted string to temporary file

        returns handle to temp file

    '''

    tmp_fh = tempfile.NamedTemporaryFile(suffix='.fa', delete=False)
    print >>tmp_fh, fasta_str
    print >>sys.stderr, 'Created temporary fasta "%s"' % tmp_fh.name
    tmp_fh.close()

    return tmp_fh

def annotate_positions(seqrec):
    ''' add position numbering to seqrec

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

def pop_row(aln, seqid):
    ''' pop a row from an alignment by sequence id

        from Bio import SeqIO
    '''

    aln_d = SeqIO.to_dict(aln)
    seq = aln_d[seqid]
    del aln_d[seqid]
    aln = Align.MultipleSeqAlignment(aln_d.itervalues())

    return seq, aln

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
    chain_ref, aln_ref = pop_row(aln_ref, chain_ref_id)
    chain_ref = annotate_positions(chain_ref)

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
    model = open_model(options['pdb_file'])
    chain = model[options['chain_id']]
    chain_seqr = chain_to_SeqRecord(chain, options['chain_id'])
    chain_seqr_m = SeqRecord.SeqRecord(chain_seqr.seq, id = mangle_id(chain_seqr.id),
                                       name = '', description = '')

    if 'refid' in options:
        # align chain to reference in alignment
        aln_l = [ seqr for seqr in aln if seqr.id == options['refid'] ]
        if len(aln_l) < 1:
            print >>sys.stderr, "refid <%s> not in reference alignment %s" % \
                                    (options['refid'], options['aln.fa'])
            sys.exit(1)
        ref_seqr = aln_l[0]
        if 'ex_aln' in options:
            ex_aligner = make_ex_aligner(PAIR_ALIGNER_CMD)
        else:
            ex_aligner = needle_align
        aln_ref = pair_align(chain_seqr_m, ref_seqr, ex_aligner)
        col_resn_aa = col_to_resn(Align.MultipleSeqAlignment([ref_seqr]),
                                  chain_seqr, aln_ref, chain_seqr_m.id)
    else:
        # align chain to alignment
        if 'ex_aln' in options:
            ex_aligner = make_ex_aligner(PROFILE_ALIGNER_CMD)
        else:
            ex_aligner = muscle_profile_align
        aln_ref = profile_align(chain_seqr_m, options['aln.fa'], ex_aligner)
        col_resn_aa = col_to_resn(aln, chain_seqr, aln_ref, chain_seqr_m.id)

    for c_r_a in col_resn_aa:
        print '\t'.join(map(str, c_r_a))

