# Sequence formatting

## Convert from fasta to phylip and PSICOV

```python
%%bash

fasta_to_phy.py < vif.fa > vif.phy
fasta_to_psicov.py < vif.fa > vif.psi
```

## Split and join alignments

```python
%%bash

join_fastas.py vif.fa A3G.fa > vif_A3G.fa
split_faa_on_col.py vif_A3G.fa 324 left.fa right.fa
```

```python
%%bash

infoalign vif.fa -auto -stdout | head
```

```python
%%bash

infoalign vif_A3G.fa -auto -stdout | head
```

```python
%%bash

infoalign left.fa -auto -stdout | head
```
