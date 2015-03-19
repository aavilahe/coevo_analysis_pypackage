# fmt_utils
An assortment of handy scripts used to format
alignments prior to a coevolution analysis

## Usage

### `concatenate_fastas.py`
Concatenates pairs of sequences from two alignments.
Headers must exactly match (seq id and comment) and be unique within each fasta.
Unpaired sequences are skipped, whitespace is removed from sequence.

```bash
usage: ./concatenate_fastas.py left.fa right.fa > left_right.fa
```

### `fasta_to_phy.py`
Convert fasta to "strict" phylip:

- header: `' Nseqs Ncol'`
- one sequence per line
	- first 8 chars are truncated seq ids
	- chars 9,10 are spaces
	- 11th char is first alignment column

Fasta header line is split on whitespace, comment is discarded.
Whitespace is removed from sequence.

```bash
usage: ./fasta_to_phy.py < fasta > phy
```

### `fasta_to_psicov.py`
Convert fasta to PSICOV readable format
	- one sequence per line, no headers
Discard all headers, remove whitespace from sequence.

```bash
usage: ./fasta_to_psicov.py < fasta > psicov
```

### `split_faa_on_col.py`
Halve an aligned fasta on a specified column. Does not modify sequence order. Appends `'_L|R'` to
sequence ids, comment is included in output untouched.

```bash
# col: 1,2,3,...,N,N+1,...,K
#       <--Left--| |--Right-->

usage: ./split_faa_on_col.py ali.fa N left.fa right.fa
```

## Author ##
Aram Avila-Herrera (Aram.Avila-Herrera@ucsf.edu)

