#!/usr/bin/env python
''' min_dists.py -- read two `get_dists.py' outputs, get min dists per res-pair

'''

import sys

def read_dist(fh):
	''' read output from get_dists.py

		left_num, right_num, left_aa, right_aa, dist
	
	'''

	dists = dict()
	for line in fh:
		left_num, right_num, left_aa, right_aa, dist = line.rstrip().split()[:5]
		dists[(int(left_num), int(right_num))] = (left_aa, right_aa, float(dist))
	
	return dists

def min_record(rec1, rec2):
	''' compare distances in record 1 and 2, return record with min distance

		rec := (left_aa, right_aa, dist)

	'''

	if rec1[2] < rec2[2]:
		return rec1
	return rec2

def print_min(dists1, dists2):
	''' print min dists for common res-pairs

	'''

	all_res_pairs = sorted(set(dists1.iterkeys()) | set(dists2.iterkeys()))

	for res_pair in all_res_pairs:
		if res_pair in dists1 and res_pair in dists2:
			min_rec = min_record(dists1[res_pair], dists2[res_pair])
		elif res_pair in dists1 and res_pair not in dists2:
			min_rec = dists1[res_pair]
		else:
			min_rec = dists2[res_pair]
		print '%d\t%d\t%s\t%s\t%.6f' % (res_pair + min_rec)

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print >>sys.stderr, 'usage: %s dists1 dists2 > mindists'
		sys.exit(1)
	
	dfn1 = open(sys.argv[1])
	dfn2 = open(sys.argv[2])

	dists1 = read_dist(dfn1)
	dists2 = read_dist(dfn2)

	print_min(dists1, dists2)

