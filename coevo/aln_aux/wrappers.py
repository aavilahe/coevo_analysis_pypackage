#!/usr/bin/python env
''' wrappers.py -- miscellaneous wrappers for dealing with alignments

'''

from os import remove
from subprocess import Popen, PIPE
from StringIO import StringIO

from Bio import Align, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Emboss.Applications import NeedleCommandline
from Bio.pairwise2 import align
from Bio.SubsMat import MatrixInfo

from .aux import make_tmp_fa, ungap_SeqRecord


def make_external_aligner(EX_ALN_CMD):
    ''' Returns a function to align to fasta files with an external aligner.

        EX_ALN_CMD: command line containing two '%s' formatting operators

        For example:
                    'muscle -profile -in1 %s -in2 %s -out /dev/stdout'
                    'needle -auto -asequence %s -bsequence %s -stdout -aformat3 fasta'

    '''

    def ex_aligner(fa1, fa2):
        ''' Aligns sequences in two fasta files with an external aligner

            fa1, fa2: filenames of fasta files to align

            Returns a MultipleSeqAlignment object

        '''

        cmd = (EX_ALN_CMD % (fa1, fa2)).split()
        cmd_sout = Popen(cmd, stdout=PIPE).communicate()[0]  # watch out for large alignments
        exaln = AlignIO.read(StringIO(cmd_sout), format = "fasta")

        return exaln

    return ex_aligner

def muscle_profile_align(fa1, fa2):
    ''' Uses muscle to profile-align two fastas

        fa1, fa2: filenames of fastas to profile-align. Must exist on disk when
                  command is called

        Returns a MultipleSeqAlignment object

    '''

    muscle_cmd = MuscleCommandline(in1 = fa1, in2 = fa2,
                                   profile = True
                                   )
    exaln = AlignIO.read(StringIO(muscle_cmd()[0]), format = "fasta")

    return exaln

def needle_align(fa1, fa2, gapopen = 10.0, gapextend = 0.5):
    ''' Uses needle to align two fastas with default gap penalties

        fa1, fa2: filenames of fastas to pairwise align. Must exist on disk when
                  command is called.
        gapopen: gap open penalty [default = 10.0]
        gapextend: gap extend penalty [default = 0.5]

        Returns a MultipleSeqAlignment object

    '''

    needle_cmd = NeedleCommandline(asequence = fa1, bsequence = fa2,
                                   outfile='/dev/stdout', aformat = 'fasta',
                                   gapopen = gapopen, gapextend = gapextend
                                   )
    exaln = AlignIO.read(StringIO(needle_cmd()[0]), format = 'fasta')

    return exaln

def profile_align_SeqRecord_to_fa(seqr, aln_fn, ex_aligner = muscle_profile_align):
    ''' Profile align SeqRecord to alignment in fasta file using external aligner.

        seqr: SeqRecord to profile-align
        aln_fn: filename of fasta to profile-align
        ex_aligner: helper function that profile-aligns sequences in two files
        
        Returns a MultipleSeqAlignment object

    '''

    tmp_fa = make_tmp_fa(seqr)
    exaln = ex_aligner(tmp_fa.name, aln_fn)
    remove(tmp_fa.name)

    return exaln

def pair_align_SeqRecords(seqr_a, seqr_b, ex_aligner = needle_align):
    ''' Pairwise align two SeqRecords using external or internal aligner.

        seqr_a, seqr_b: SeqRecords to align
        ex_aligner: helper function that aligns sequences in two files
                    *If ex_aligner is None, use internal aligner*

        Internal aligner:
            Bio.pairwise2.align.globalds() with Bio.SubsMat.MatrixInfo.Blosum62
            and default gap penalties (gapopen = -10.0, gapextend = -0.5)

        Returns a MultipleSeqAlignment object

    '''

    if ex_aligner is None:
        inaln = align.globalds(ungap_SeqRecord(seqr_a), ungap_SeqRecord(seqr_b),
                               MatrixInfo.blosum62, -10.0, -0.5)
        exaln = Align.MultipleSeqAlignment(inaln[0][:2])
    else:
        tmp_fa = make_tmp_fa(ungap_SeqRecord(seqr_a))
        tmp_ref_fa = make_tmp_fa(ungap_SeqRecord(seqr_b))
        exaln = ex_aligner(tmp_fa.name, tmp_ref_fa.name)
        remove(tmp_fa.name)
        remove(tmp_ref_fa.name)

    return exaln


