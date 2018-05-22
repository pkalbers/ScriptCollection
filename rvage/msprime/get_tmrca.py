import csv
import sys
import msprime as msp


name =  sys.argv[1]
hdf5 =  sys.argv[2]

data = msp.load(hdf5)

input = open(name, 'r')
lines = csv.DictReader(input, delimiter=" ")


print("MarkerID Clock Fk SampleID0 Chr0 SampleID1 Chr1 Shared SegmentLHS SegmentRHS Pass Shape Rate MRCA TMRCA")


def run_group(idx, grp, data):
	for tree in data.trees():
		for mut in tree.mutations():
			if idx > mut.index:
				continue
			if idx == mut.index:
				for item in grp:
					mid = int(item['MarkerID'])
					mfk = int(item['Fk'])
					shr = int(item['Shared'])
					id0 = (int(item['SampleID0']) * 2) + int(item['Chr0'])
					id1 = (int(item['SampleID1']) * 2) + int(item['Chr1'])
					mrca = tree.get_mrca(id0, id1)
					time = tree.get_time(mrca)
					print("%d %s %d %s %s %s %s %d %s %s %s %s %s %d %.6f" % (mid, item['Clock'], mfk, item['SampleID0'], item['Chr0'], item['SampleID1'], item['Chr1'], shr, item['SegmentLHS'], item['SegmentRHS'], item['Pass'], item['Shape'], item['Rate'], mrca, time))
				return


idx = -1
grp = []

for line in lines:
	jdx = int(line['MarkerID'])
	if idx != jdx:
		if len(grp) > 0:
			run_group(idx, grp, data)
		grp = []
		idx = jdx
	grp.append(line)

if len(grp) > 0:
	run_group(idx, grp, data)
