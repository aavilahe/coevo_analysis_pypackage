# Mapping structure to alignment and distances example

An example workflow for getting distances between
alignment columns using a PDB structure.

## Map 3DGE---An A2B2 heterocomplex to an alignment of ortholog pairs.

We will map [3DGE](http://www.rcsb.org/pdb/explore/explore.do?structureId=3DGE)
to the HisKA-RR alignment used in [Procaccini *et al* PLoS One
2011](doi.org/10.1371/journal.pone.0019729) and [Avila-Herrera & Pollard BMC
Binf 2015](http://doi.org/10.1186/s12859-015-0677-y).

##### Splitting

We must first split the concatenated alignment into its two domains on column
87.

```python
%%bash

split_faa_on_col.py HisKA_RR.fa 87 HisKA.fa RR.fa
```

##### Map alignment columns to residue numberings (resnum)

Chain A has contacts with chains C and D. And chains C and D are identical (but
we can never be too sure, so we map anyways).

*`map_column_to_resnum.py` works on gzipped PDB files, too!*

```python
%%bash

map_column_to_resnum.py A 3DGE.pdb.gz HisKA.fa > HisKA-A.tsv
map_column_to_resnum.py C 3DGE.pdb.gz RR.fa > RR-C.tsv
map_column_to_resnum.py D 3DGE.pdb.gz RR.fa > RR-D.tsv
```

```python
%%bash

head HisKA-A.tsv RR-C.tsv RR-D.tsv
```

## Get atomic distances between HisKA and RR residues

`get_dists.py` calculates distances between chains and uses the user-supplied
column-to-resnum map to report distances between alignment columns.

Remember `RR.fa` maps to both chains C and D.

```python
%%bash

get_dists.py --chainL=A --chainR=C --mapL=HisKA-A.tsv --mapR=RR-C.tsv --dist_atoms=Cb 3DGE.pdb.gz > 3DGE.A_C.Cb_dist.tsv
get_dists.py --chainL=A --chainR=D --mapL=HisKA-A.tsv --mapR=RR-D.tsv --dist_atoms=Cb 3DGE.pdb.gz > 3DGE.A_D.Cb_dist.tsv

```

```python
%%bash

head 3DGE.A_C.Cb_dist.tsv 3DGE.A_D.Cb_dist.tsv
```

##### Get minimum distance for each pair of columns

Because each column in `RR.fa` maps to two chains, we will choose the smallest
distance for each pair of HisKA-RR columns.

```python
%%bash

min_dists.py 3DGE.A_C.Cb_dist.tsv 3DGE.A_D.Cb_dist.tsv > 3DGE.Cb_dist.tsv
```

```python
%%bash

head 3DGE.Cb_dist.tsv
```

