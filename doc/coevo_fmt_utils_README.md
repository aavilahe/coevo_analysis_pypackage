# Sequence formatting

An assortment of handy scripts used to format
alignments prior to a coevolution analysis.

## `join_fastas.py`

Joins two fasta alignments on common sequence identifiers (must exactly match
and be unique within each fasta).

```bash

usage: join_fastas.py left.fa right.fa > left_right.fa

```

## `fasta_to_phy.py`

Convert fasta to "stricter" phylip:

- Header: `' Nseqs Ncol'`
- One sequence per line
	- First 8 chars are truncated seq ids
	- Characters 9 and 10 are spaces
	- 11th character is the first alignment column

Fasta header line is split on whitespace, comment is discarded.

```bash

usage: fasta_to_phy.py < aln.fa > aln.phy

```

## `fasta_to_psicov.py`

Convert fasta to PSICOV readable format
	- One sequence per line
    - Discards all headers

```bash

usage: fasta_to_psicov.py < aln.fa > aln.psi

```

## `split_faa_on_col.py`

Split an aligned fasta on a specified column.
Headers are untouched.

```bash

# col: 1,2,3,...,N,N+1,...,K
#       <--Left--| |--Right-->

usage: split_faa_on_col.py aln.fa N left.fa right.fa

```

