#!/usr/bin/env python
''' get_dists.py -- gets distances between residues in two or within same chain

	user can specify residue-residue distance type:
		- beta carbon (default)
		- nearest nonhydrogen atom
		- nearest atom
	
	depends:
		- biopython >= 1.64
		- numpy
	
	output:
		resnum_left resnum_right resname_left resname_right distance

'''

import sys
import getopt
import gzip
from itertools import product, combinations
from os.path import basename
from re import sub
from numpy.linalg import norm
import Bio.SeqUtils
import Bio.PDB

__author__ = 'Aram Avila-Herrera'

def getCBCoords(residue):
	''' get beta carbon coordinates from residue

	'''

	if residue.get_resname() == 'GLY':
		return [ residue['CA'].coord ]
	else:
		return [ residue['CB'].coord ]

def getNoHCoords(residue):
	''' get non-hydrogen atom coordinates from residue

	'''

	return [ atom.coord for atom in residue.get_list() if not atom.get_id().startswith('H') ]

def getAllCoords(residue):
	''' get all atom coordinates from residue

	'''

	return [ atom.coord for atom in residue.get_list() ]

def getNonHetResidues(chain):
	''' get list of residues from chain, remove those with het flag

	'''

	all_res = chain.get_residues()
	nonhet_res = (res for res in all_res if res.get_id()[0] == ' ')
	return nonhet_res

def hasStructureCarbon(res):
	''' returns true if residue has valid CB (CA if res is glycine)\
	
	'''

	if res.get_resname() == 'GLY':
		if 'CA' in res:
			return True
		return False
	elif 'CB' in res:
		return True
	return False

def getStructResidues(nonhet_res):
	''' keep only residues with CB coordinates (or CA if glycine)

	'''

	struct_nonhet_res_residues = ( res for res in nonhet_res if hasStructureCarbon(res) )
	return struct_nonhet_res_residues

def calcResidueDistance(residue_left, residue_right, dtype='Cb'):
	''' get euclidean distance between residues

		dtype := [ 'Cb' | 'NoH' | 'Any' ]

	'''
	
	if dtype == 'Cb':
		getCoords = getCBCoords
	elif dtype == 'NoH':
		getCoords = getNoHCoords
	else:
		getCoords = getAllCoords
	
	coords_left = getCoords(residue_left)
	coords_right = getCoords(residue_right)
	
	dists = ( norm(X1-X2) for (X1, X2) in product(coords_left, coords_right) )
	
	# interested in closest atom-atom distance between residues
	return min(dists)

def getInterChainDistances(chain_left, chain_right, dtype='Cb'):
	''' get distances for all pairs of residues between two chains
		
		returns (resnumL, resnumR, aaL, aaR, dist)

	'''

	res_pairs = product(getStructResidues(getNonHetResidues(chain_left)), getStructResidues(getNonHetResidues(chain_right)))
	return [
				(
					rl.id[1], rr.id[1],                                           # resnums
					Bio.SeqUtils.seq1(rl.resname), Bio.SeqUtils.seq1(rr.resname), # 3->1 letter resnames
					calcResidueDistance(rl, rr, dtype)                            # distance
				) for (rl, rr) in res_pairs
			]

def getIntraChainDistances(chain, dtype='Cb'):
	''' get distances for all pairs of residues within one chain

		returns (resnumL, resnumR, aaL, aaR, dist)

	'''

	res_pairs = combinations(getStructResidues(getNonHetResidues(chain)), r=2)
	return [
				(
					rl.id[1], rr.id[1],                                           # resnums
					Bio.SeqUtils.seq1(rl.resname), Bio.SeqUtils.seq1(rr.resname), # 3->1 letter resnames
					calcResidueDistance(rl, rr, dtype)                            # distance
				) for (rl, rr) in res_pairs
			]

def commandLineOpts(args):
	''' Parse command line or config file for options and return a dict.

		using getopt for compatibility

	'''

	optlist, args = getopt.getopt(args, 'hc:d:m:n:',
					[ 
						'help',
						'chainL=',
						'chainR=',
						'mapL=',
						'mapR=',
						'dtype='
					])
	options = dict()
	usage = (
						'usage: %s [ options ] pdb_file\n\n' +\
						'Options:\n'+\
						'\t-h, --help                     show this help message and exit\n' +\
						'\t-c, --chainL <chain id>        (required) specify single or left chain id\n' +\
						'\t--chainR <chain id>            specify right chain id\n' +\
						'\t-m, --mapL <left.col2resnum>   specify column to resnum mapping file for left chain\n' +\
						'\t-n, --mapR <right.col2resnum>  specify column to resnum mapping file for right chain\n' +\
						'\t-d, --dtype Cb | NoH | Any   specify beta carbon (default), non-hydrogen, or any atom distances\n'
			) % sys.argv[0]
	
	# defaults
	options['dtype'] = 'Cb'

	# assign options
	for opt, val in optlist:
		if opt in ('-h', '--help'):
			sys.exit(usage)
		if opt in ('-c', '--chainL'):
			options['chainL'] = val
		if opt == '--chainR':
			options['chainR'] = val
		if opt in ('-m', '--mapL'):
			options['mapL'] = val
		if opt in ('-n', '--mapR'):
			options['mapR'] = val
		if opt in ('-d', '--dtype'):
			options['dtype'] = val
	
	# check arguments and required options
	
	if len(args) < 1:
		err = 'Error: specify a pdb_file'
		sys.exit(err + '\n' + usage)
	options['pdb_file'] = args[0]
	
	if 'chainL' not in options:
		err = 'Error: specify a single or left chain'
		sys.exit(err + '\n' + usage)

	return options

def openModel(pdb_fn):
	''' opens first model in given pdb file name

		if file name ends in .gz, will attempt to open as a gzipped file

	'''

	if pdb_fn.endswith('.gz'):
		pdb_fh = gzip.open(pdb_fn)
	else:
		pdb_fh = open(pdb_fn)

	pdb_id = basename(pdb_fn).split('.')[0]
	structure = Bio.PDB.PDBParser().get_structure(pdb_id, pdb_fh)
	model = structure[0] # assume 1 model in pdb file
	
	return model

def readMap(fn):
	''' read map from resnum to column

	'''

	fh = open(fn)
	rn2col = dict(( map(int, reversed(line.split()[:2])) for line in fh.readlines() ))

	return rn2col

def convertIdx(results_list, mapdict, chain_sel):
	''' converts residue numberings using mapdict

		results_list := list of tuples (resnum_L, resnum_R, aa_L, aa_R, distance)
		mapdict := dict of mapping eg. {key = resnum, value = colnum}
		chain_sel := int specifying which item in tuple to modify (eg. 0 = resnum_L, 1 = renum_R)
	
	'''

	# mapped_res = [ tup[:chain_sel] + (mapdict[tup[chain_sel]],) + tup[chain_sel+1:] for tup in results_list if tup[chain_sel] in mapdict ]

	mapped_res = list()
	for (i, tup) in enumerate(results_list):
		idx_from = tup[chain_sel]
		if idx_from in mapdict:
			idx_to = mapdict[idx_from]
			newtup = tup[:chain_sel] + tuple([idx_to]) + tup[chain_sel+1:]
			mapped_res.append(newtup)
	
	return mapped_res
	
def convertResnums(results, options):
	''' if mapping files given, convert resnums to alignment columns
		
	'''

	if 'mapL' in options:
		mapL = readMap(options['mapL'])
		results = convertIdx(results, mapL, 0)
	if 'mapR' in options:
		mapR = readMap(options['mapR'])
		results = convertIdx(results, mapR, 1)

	return results

if __name__ == "__main__":
	options = commandLineOpts(sys.argv[1:])
	model = openModel(options['pdb_file'])

	if 'chainR' in options:
		results = getInterChainDistances(model[options['chainL']], model[options['chainR']])
	else:
		results = getIntraChainDistances(model[options['chainL']])

	#results = [ map(str, result) for result in results ]
	results = convertResnums(results, options)
	for result in results:
		#print '\t'.join(result)
		print '%s\t%s\t%s\t%s\t%.6f' % tuple(result)
	
