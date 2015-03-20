#!/bin/bash

../fasta_to_phy.py < vif.fa > vif.phy
../fasta_to_psicov.py < vif.fa > vif.psi
../join_fastas.py vif.fa A3G.fa > vif_A3G.fa
../split_faa_on_col.py vif_A3G.fa 324 left.fa right.fa
