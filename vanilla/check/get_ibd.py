import csv
import sys
import msprime as msp

name =  sys.argv[1]
hdf5 =  sys.argv[2]

file = open(name, 'r')
list = csv.DictReader(file, delimiter=" ")

data = msp.load(hdf5)

print("MarkerID SampleID0 Chr0 SampleID1 Chr1 Shared TrueLHS TrueRHS")

for line in list:
	mid = int(line['MarkerID'])
	id0 = int(line['SampleID0']) * 2
	id1 = int(line['SampleID1']) * 2
	shr = int(line['Shared'])

	if line['Chr0'] == 'P':
		id0 = id0 + 1
	if line['Chr1'] == 'P':
		id1 = id1 + 1

	flag = False
	ends = False
	last = -1
	xlhs = -1
	xrhs = -1

	cx = 0

	for tree in data.trees():
		if ends:
			break
		for mut in tree.mutations():
			mrca = tree.get_mrca(id0, id1)

			if flag:
				xrhs = mut.index

			if last != mrca:
				if flag:
					ends = True
					break
				else:
					xlhs = cx

			if mid == mut.index:
				flag = True
				xrhs = mut.index

			last = mrca
			cx = mut.index

	print("%d %s %s %s %s %d %d %d" % (mid, line['SampleID0'], line['Chr0'], line['SampleID1'], line['Chr1'], shr, xlhs, xrhs))
