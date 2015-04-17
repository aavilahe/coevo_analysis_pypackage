#!/bin/bash
# Example workflow for getting distances between alignment columns
# using structural models
#
#
# 3DGE - A2B2 heterocomplex

BIN="../../bin"

PDB='3DGE.pdb.gz'

CHAIN_L1='A'
#CHAIN_L2='B' #ignored for now
CHAIN_R1='C'
CHAIN_R2='D'

FA_L='HisKA.fa'
FA_R='RR.fa'

# map columns to alignment (note chain A has contacts with C and D) 
MAPL1="${FA_L%.fa}-${CHAIN_L1}.col_res_aa"
MAPR1="${FA_R%.fa}-${CHAIN_R1}.col_res_aa"
MAPR2="${FA_R%.fa}-${CHAIN_R2}.col_res_aa"

python ${BIN}/map_column_to_resnum.py --user_aln ${CHAIN_L1} ${PDB} ${FA_L} > ${MAPL1}
python ${BIN}/map_column_to_resnum.py --user_aln ${CHAIN_R1} ${PDB} ${FA_R} > ${MAPR1}
python ${BIN}/map_column_to_resnum.py --user_aln ${CHAIN_R2} ${PDB} ${FA_R} > ${MAPR2}

# get resn-resn distances
#RDIST1="${PDB%.pdb*}.${CHAIN_L1}_${CHAIN_R1}.Cb.res_dist"
#RDIST2="${PDB%.pdb*}.${CHAIN_L1}_${CHAIN_R2}.Cb.res_dist"

#python ${BIN}/get_dists.py --chainL=${CHAIN_L1} --chainR=${CHAIN_R1} --dist_atoms=Cb ${PDB} > ${RDIST1}
#python ${BIN}/get_dists.py --chainL=${CHAIN_L1} --chainR=${CHAIN_R2} --dist_atoms=Cb ${PDB} > ${RDIST2}

# get col-col distances
CDIST1="${PDB%.pdb*}.${CHAIN_L1}_${CHAIN_R1}.Cb.col_dist"
CDIST2="${PDB%.pdb*}.${CHAIN_L1}_${CHAIN_R2}.Cb.col_dist"

python ${BIN}/get_dists.py --chainL=${CHAIN_L1} --chainR=${CHAIN_R1} --mapL=${MAPL1} --mapR=${MAPR1} --dist_atoms=Cb ${PDB} > ${CDIST1}
python ${BIN}/get_dists.py --chainL=${CHAIN_L1} --chainR=${CHAIN_R2} --mapL=${MAPL1} --mapR=${MAPR2} --dist_atoms=Cb ${PDB} > ${CDIST2}

# convert resn to column
#python ${BIN}/convert_resnums_to_columns.py ${RDIST1} ${MAPL1} ${MAPR1} > ${CDIST1}
#python ${BIN}/convert_resnums_to_columns.py ${RDIST2} ${MAPL1} ${MAPR2} > ${CDIST2}

# get min distance
MDIST="${PDB%.pdb*}.Cb.mindist"
python ${BIN}/min_dists.py ${CDIST1} ${CDIST2} > ${MDIST}
