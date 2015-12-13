# Sequence formatting

## Convert from fasta to phylip and PSICOV

```{.python .input  n=1}
%%bash

fasta_to_phy.py < vif.fa > vif.phy
fasta_to_psicov.py < vif.fa > vif.psi
```

## Split and join alignments

```{.python .input  n=2}
%%bash

join_fastas.py vif.fa A3G.fa > vif_A3G.fa
split_faa_on_col.py vif_A3G.fa 324 left.fa right.fa
```

```{.python .input  n=3}
%%bash

infoalign vif.fa -auto -stdout | head
```

```
# USA             Name        SeqLen	AlignLen	Gaps	GapLen	Ident	Similar	Differ	% Change	Weight	Description
fasta::vif.fa:0_HIV1_h	0_HIV1_h      192	324	20	132	75	19	98	76.851852	1.000000	
fasta::vif.fa:1_HIV2_h	1_HIV2_h      215	324	17	109	96	26	93	70.370369	1.000000	
fasta::vif.fa:2_SIVcpz	2_SIVcpz      192	324	20	132	79	22	91	75.617287	1.000000	
fasta::vif.fa:3_SIVmac	3_SIVmac      214	324	17	110	106	21	87	67.283951	1.000000	
fasta::vif.fa:4_SIVcol	4_SIVcol      159	324	20	165	37	18	104	88.580246	1.000000	
fasta::vif.fa:5_SIVdeb	5_SIVdeb      250	324	15	74	70	34	146	78.395065	1.000000	
fasta::vif.fa:6_SIVgor	6_SIVgor      192	324	20	132	75	22	95	76.851852	1.000000	
fasta::vif.fa:7_SIVgri	7_SIVgri      219	324	16	105	85	30	104	73.765434	1.000000	
fasta::vif.fa:8_SIVmne	8_SIVmne      214	324	17	110	105	22	87	67.592590	1.000000	

```

```{.python .input  n=4}
%%bash

infoalign vif_A3G.fa -auto -stdout | head
```

```
# USA             Name        SeqLen	AlignLen	Gaps	GapLen	Ident	Similar	Differ	% Change	Weight	Description
fasta::vif_A3G.fa:0_HIV1_h	0_HIV1_h      576	712	22	136	384	42	150	46.067417	1.000000	
fasta::vif_A3G.fa:1_HIV2_h	1_HIV2_h      599	712	19	113	405	49	145	43.117977	1.000000	
fasta::vif_A3G.fa:2_SIVcpz	2_SIVcpz      576	712	22	136	388	47	141	45.505619	1.000000	
fasta::vif_A3G.fa:3_SIVmac	3_SIVmac      598	712	19	114	456	32	110	35.955055	1.000000	
fasta::vif_A3G.fa:4_SIVcol	4_SIVcol      544	712	23	168	370	31	143	48.033707	1.000000	
fasta::vif_A3G.fa:5_SIVdeb	5_SIVdeb      633	712	18	79	433	43	157	39.185394	1.000000	
fasta::vif_A3G.fa:6_SIVgor	6_SIVgor      576	712	22	136	384	48	144	46.067417	1.000000	
fasta::vif_A3G.fa:7_SIVgri	7_SIVgri      602	712	19	110	452	39	111	36.516853	1.000000	
fasta::vif_A3G.fa:8_SIVmne	8_SIVmne      590	712	21	122	449	34	107	36.938202	1.000000	

```

```{.python .input  n=5}
%%bash

infoalign left.fa -auto -stdout | head
```

```
# USA             Name        SeqLen	AlignLen	Gaps	GapLen	Ident	Similar	Differ	% Change	Weight	Description
fasta::left.fa:0_HIV1_h	0_HIV1_h      192	324	20	132	75	19	98	76.851852	1.000000	
fasta::left.fa:1_HIV2_h	1_HIV2_h      215	324	17	109	96	26	93	70.370369	1.000000	
fasta::left.fa:2_SIVcpz	2_SIVcpz      192	324	20	132	79	22	91	75.617287	1.000000	
fasta::left.fa:3_SIVmac	3_SIVmac      214	324	17	110	106	21	87	67.283951	1.000000	
fasta::left.fa:4_SIVcol	4_SIVcol      159	324	20	165	37	18	104	88.580246	1.000000	
fasta::left.fa:5_SIVdeb	5_SIVdeb      250	324	15	74	70	34	146	78.395065	1.000000	
fasta::left.fa:6_SIVgor	6_SIVgor      192	324	20	132	75	22	95	76.851852	1.000000	
fasta::left.fa:7_SIVgri	7_SIVgri      219	324	16	105	85	30	104	73.765434	1.000000	
fasta::left.fa:8_SIVmne	8_SIVmne      214	324	17	110	105	22	87	67.592590	1.000000	

```
