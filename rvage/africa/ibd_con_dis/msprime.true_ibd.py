import csv
import sys
import msprime as msp

name =  sys.argv[1]
hdf5 =  sys.argv[2]

file = open(name, 'r')
list = csv.DictReader(file, delimiter=" ")

data = msp.load(hdf5)

print("MarkerID Fk SampleID0 Chr0 Act0 SampleID1 Chr1 Act1 Shared LHS RHS TrueLHS TrueRHS MRCA TMRCA")

runs = []
count = 0
for line in list:
	line['mid'] = int(line['MarkerID'])
	line['mfk'] = int(line['Fk'])
	line['id0'] = int(line['SampleID0']) * 2
	line['id1'] = int(line['SampleID1']) * 2
	line['shr'] = int(line['Shared'])
	line['infL'] = int(line['SegmentLHS'])
	line['infR'] = int(line['SegmentRHS'])

	if line['Act0'] == '1':
		line['id0'] = line['id0'] + 1
	if line['Act1'] == '1':
		line['id1'] = line['id1'] + 1

	line['flag'] = False
	line['ends'] = False
	line['last'] = -1
	line['time'] = -1.0
	line['xlhs'] = -1
	line['xrhs'] = -1

	line['cx'] = 0

	count = count + 1

	runs.append(line)

	if count == 50:
		for tree in data.trees():
			a = 0
			b = 0
			for run in runs:
				a = a + 1
				if run['ends']:
					b = b + 1
			if b > 0 and a == b:
				break
			for mut in tree.mutations():
				for run in runs:
					if run['ends']:
						continue

					mrca = tree.get_mrca(run['id0'], run['id1'])
					tmrca = tree.get_time(mrca)

					if run['flag']:
						run['xrhs'] = mut.index

					if run['last'] != mrca:
						if run['flag']:
							run['ends'] = True
						else:
							run['xlhs'] = run['cx']

					if run['mid'] == mut.index:
						run['flag'] = True
						run['xrhs'] = mut.index

					run['last'] = mrca
					run['time'] = tmrca
					run['cx'] = mut.index

		for run in runs:
			print("%d %d %s %s %s %s %s %s %d %d %d %d %d %d %.6f" % (run['mid'], run['mfk'], run['SampleID0'], run['Chr0'], run['Act0'], run['SampleID1'], run['Chr1'], run['Act1'], run['shr'], run['infL'], run['infR'], run['xlhs'], run['xrhs'], run['last'], run['time']))

		runs = []
		count = 0

