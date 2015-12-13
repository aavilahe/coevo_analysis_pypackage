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

```{.python .input  n=1}
%%bash

split_faa_on_col.py HisKA_RR.fa 87 HisKA.fa RR.fa
```

##### Map alignment columns to residue numberings (resnum)

Chain A has contacts with chains C and D. And chains C and D are identical (but
we can never be too sure, so we map anyways).

*`map_column_to_resnum.py` works on gzipped PDB files, too!*

```{.python .input  n=2}
%%bash

map_column_to_resnum.py A 3DGE.pdb.gz HisKA.fa > HisKA-A.tsv
map_column_to_resnum.py C 3DGE.pdb.gz RR.fa > RR-C.tsv
map_column_to_resnum.py D 3DGE.pdb.gz RR.fa > RR-D.tsv
```

```
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6530.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6535.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6540.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6567.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6607.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6644.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6704.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6731.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain N is discontinuous at line 6767.
  PDBConstructionWarning)
Created temporary fasta "/var/folders/3p/5y9jnynj7ns0m_hs_zsgd5lc0000gn/T/tmpq_LSzY.fa"
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6530.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6535.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6540.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6567.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6607.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6644.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6704.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6731.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain N is discontinuous at line 6767.
  PDBConstructionWarning)
Created temporary fasta "/var/folders/3p/5y9jnynj7ns0m_hs_zsgd5lc0000gn/T/tmpFlf_x7.fa"
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6530.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6535.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6540.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6567.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6607.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6644.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6704.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6731.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain N is discontinuous at line 6767.
  PDBConstructionWarning)
Created temporary fasta "/var/folders/3p/5y9jnynj7ns0m_hs_zsgd5lc0000gn/T/tmphJB4RY.fa"

```

```{.python .input  n=3}
%%bash

head HisKA-A.tsv RR-C.tsv RR-D.tsv
```

```
==> HisKA-A.tsv <==
Column	resn	AA
0	250	M
1	251	K
2	252	T
3	253	E
4	254	F
5	255	I
6	256	A
7	257	N
8	258	I

==> RR-C.tsv <==
Column	resn	AA
0	3	K
1	4	K
2	5	V
3	6	L
4	7	L
5	8	V
6	9	D
7	10	D
8	11	S

==> RR-D.tsv <==
Column	resn	AA
0	3	K
1	4	K
2	5	V
3	6	L
4	7	L
5	8	V
6	9	D
7	10	D
8	11	S

```

## Get atomic distances between HisKA and RR residues

`get_dists.py` calculates distances between chains and uses the user-supplied
column-to-resnum map to report distances between alignment columns.

Remember `RR.fa` maps to both chains C and D.

```{.python .input  n=4}
%%bash

get_dists.py --chainL=A --chainR=C --mapL=HisKA-A.tsv --mapR=RR-C.tsv --dist_atoms=Cb 3DGE.pdb.gz > 3DGE.A_C.Cb_dist.tsv
get_dists.py --chainL=A --chainR=D --mapL=HisKA-A.tsv --mapR=RR-D.tsv --dist_atoms=Cb 3DGE.pdb.gz > 3DGE.A_D.Cb_dist.tsv

```

```
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6530.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6535.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6540.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6567.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6607.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6644.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6704.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6731.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain N is discontinuous at line 6767.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6530.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6535.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6540.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6567.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 6607.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 6644.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 6704.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain D is discontinuous at line 6731.
  PDBConstructionWarning)
/usr/local/lib/python2.7/site-packages/Bio/PDB/StructureBuilder.py:86: PDBConstructionWarning: WARNING: Chain N is discontinuous at line 6767.
  PDBConstructionWarning)

```

```{.python .input  n=5}
%%bash

head 3DGE.A_C.Cb_dist.tsv 3DGE.A_D.Cb_dist.tsv
```

```
==> 3DGE.A_C.Cb_dist.tsv <==
Left_Column	Right_Column	Distance	Left_AA	Right_AA
0	0	46.907352	M	K
1	0	46.807739	K	K
2	0	44.143784	T	K
3	0	40.972687	E	K
4	0	42.922810	F	K
5	0	40.892334	I	K
6	0	36.651306	A	K
7	0	37.672855	N	K
8	0	38.631340	I	K

==> 3DGE.A_D.Cb_dist.tsv <==
Left_Column	Right_Column	Distance	Left_AA	Right_AA
0	0	53.452705	M	K
1	0	49.687504	K	K
2	0	53.047859	T	K
3	0	52.783432	E	K
4	0	48.178757	F	K
5	0	47.743767	I	K
6	0	49.812580	A	K
7	0	47.396492	N	K
8	0	43.025303	I	K

```

##### Get minimum distance for each pair of columns

Because each column in `RR.fa` maps to two chains, we will choose the smallest
distance for each pair of HisKA-RR columns.

```{.python .input  n=6}
%%bash

min_dists.py 3DGE.A_C.Cb_dist.tsv 3DGE.A_D.Cb_dist.tsv > 3DGE.Cb_dist.tsv
```

```{.python .input  n=7}
%%bash

head 3DGE.Cb_dist.tsv
```

```
Left_Column	Right_Column	Distance	Left_AA	Right_AA
0	0	46.907352	M	K
0	1	45.309952	M	K
0	2	40.697777	M	V
0	3	38.744755	M	L
0	4	36.125793	M	L
0	5	32.445393	M	V
0	6	29.931234	M	D
0	7	26.683163	M	D
0	8	27.716152	M	S

```
