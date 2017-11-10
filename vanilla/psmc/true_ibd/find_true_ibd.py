import os
import csv
import sys
import string
import random
import msprime as msp
from multiprocessing import Pool
from shutil import copyfile

name =  sys.argv[1]
hdf5 =  sys.argv[2]

core = 0
if len(sys.argv) > 3:
	core =  int(sys.argv[3])

buff = 100
if len(sys.argv) > 4:
	buff =  int(sys.argv[4])



def find_segment(line, data):
	mid = int(line['MarkerID'])
	mfk = int(line['Fk'])
	id0 = int(line['SampleID0']) * 2
	id1 = int(line['SampleID1']) * 2
	shr = int(line['Shared'])

	if line['Chr0'] == '1':
		id0 = id0 + 1

	if line['Chr1'] == '1':
		id1 = id1 + 1

	left = True
	done = False

	left_brkp = 1
	rght_brkp = -1
	prev_mrca = -1
	prev_time = -1.0
	this_mrca = -1
	this_time = -1.0

	seq = [None] * (mid + 1)

	for tree in data.trees():
		if done:
			break

		if tree.get_num_mutations() == 0:
			continue

		for mut in tree.mutations():
			mrca = tree.get_mrca(id0, id1)
			time = tree.get_time(mrca)

			if left:
				if mid < mut.index:
					this_mrca = prev_mrca
					this_time = prev_time
					left = False
				else:
					seq[mut.index] = (mrca, time)

			if not left:
				if prev_mrca != mrca or abs(prev_time - time) > 0.1:
					rght_brkp = mut.index
					done = True
					break

			prev_mrca = mrca
			prev_time = time
			break

	if left:
		this_mrca = prev_mrca
		this_time = prev_time

	if not done:
		rght_brkp = data.get_num_mutations() - 1

	last = mid
	done = False

	for i in reversed(range(mid)):
		if seq[i] == None:
			continue
		if this_mrca != seq[i][0] or abs(this_time - seq[i][1]) > 0.1:
			left_brkp = last - 1
			done = True
			break
		else:
			last = i

	if not done:
		left_brkp = 0

	out = "%d %d %s %s %s %s %d %d %d %d %.6f" % (mid, mfk, line['SampleID0'], line['Chr0'], line['SampleID1'], line['Chr1'], shr, left_brkp, rght_brkp, this_mrca, this_time)
	return out


def tw_find_segment(chunk):
	tmpf = "tmp." + ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(32)) + ".hdf5"
	copyfile(chunk['Data'], tmpf)
	data = msp.load(chunk['Data'])
	rslt = []
	for line in chunk['List']:
		rslt.append(find_segment(line, data))
	os.remove(tmpf)
	return rslt


def make_chunks(ls, n, f):
    for i in range(0, len(ls), n):
        yield { 'Data': f, 'List': ls[i:i + n] }



input = open(name, 'r')
lines = csv.DictReader(input, delimiter=" ")


items = []
for line in lines:
	items.append(line)


if core < 2:
	data = msp.load(hdf5)
	print("MarkerID Fk SampleID0 Chr0 SampleID1 Chr1 Shared LHS RHS MRCA TMRCA")
	for item in items:
		print(find_segment(item, data))
else:
	chunks = []
	for chunk in make_chunks(items, buff, hdf5):
		chunks.append(chunk)

	pool = Pool(processes=core)
	rslt = pool.map(tw_find_segment, chunks)
	pool.close()
	pool.join()

	print("MarkerID Fk SampleID0 Chr0 SampleID1 Chr1 Shared LHS RHS MRCA TMRCA")

	for rsl in rslt:
		for rs in rsl:
			print(rs)


