# Sequence to structure

These are scripts to help visualize alignment column
scores and labels onto structures (PDB files).

## `map_column_to_resnum.py`

Use this to generate a map of alignment columns to residue numbers (resnums) in
a corresponding chain in a structural model. You can use your alignment or a
sequence in your alignment as a reference.

```bash

map_column_to_resnum.py chain_id model.pdb aln.fa > col_resn_aa.map

```

The default profile aligner (for mapping your structure to your alignment) is
`muscle`.  The default pairwise global aligner is `needle` from Emboss.

If you are brave, you can edit the source in `map_column_to_resnum.py` to
specify a command that outputs an aligned fasta.

## `make_attributes.py`

Loads map of alignment columns to residue numbers (resnums) and coevolution scores in a tabular
format.
Writes chimera readable attributes file, one attribute per score column.

```

python make_attributes.py --LeftMap=ProtX_col_resnum_aa.map Scores.tsv > ProtX_attributes.txt

```
