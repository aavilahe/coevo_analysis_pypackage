# coevo_struct_viz

These are scripts to help visualize alignment column
scores and labels onto structures (PDB files).


## Usage

### `map_column_to_resnum.py`

Use this to map alignment columns to resnums in a corresponding chain
in a structural model. You can use the alignment or a sequence in the alignment as a reference.

```bash
$ ./map_column_to_resnum.py chain_id model.pdb aln.fa > col_resn_aa.tsv

```

Default profile aligner (for mapping to alignment) is `muscle`. Default pairwise global aligner
is `needle` from Emboss.

# map co-evolution scores to residues
python make_attributes.py --LeftMap=ProtX_resnum_col.map Results.tab > ProtX_attributes.txt

### `map_column_to_resnum.py` ###
#### *Mapping not optimal* ####
*Mapping through reference sequence position removes insertions
relative to reference from analysis, may not be desired* 

Takes phylip alignment as input and generates
3 mappings

1. Column to reference sequence
2. Reference sequence position to PDB chain sequence position
3. PDB chain sequence position to PDB chain resnum

And prints a map of Column to PDB chain resnum

### `make_attributes.py` ###

Loads map of Column to PDB chain resnum and tabular co-evolution results file.
Writes chimera readable attributes file, one attribute per column (skips Alignment_Column columns).




