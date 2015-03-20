# fmt_utils

An assortment of handy scripts used to format
alignments prior to a coevolution analysis.

## Usage



### `join_fastas.py`

Joins two fasta alignments on common sequence identifiers (must exactly match
and be unique within each fasta).

```bash
usage: ./join_fastas.py left.fa right.fa > left_right.fa
```

### `fasta_to_phy.py`

Convert fasta to "stricter" phylip:

- header: `' Nseqs Ncol'`
- one sequence per line
	- first 8 chars are truncated seq ids
	- chars 9,10 are spaces
	- 11th char is first alignment column

Fasta header line is split on whitespace, comment is discarded.

```bash
usage: ./fasta_to_phy.py < aln.fa > aln.phy
```

### `fasta_to_psicov.py`

Convert fasta to PSICOV readable format
	- one sequence per line
    - discards all headers

```bash
usage: ./fasta_to_psicov.py < aln.fa > aln.psi
```

### `split_faa_on_col.py`
Split an aligned fasta on a specified column.
Headers are untouched.

```bash
# col: 1,2,3,...,N,N+1,...,K
#       <--Left--| |--Right-->

usage: ./split_faa_on_col.py aln.fa N left.fa right.fa
```

## Depends

- Biopython

## Author

Aram Avila-Herrera (Aram.Avila-Herrera@ucsf.edu)

