import csv
import sys
import msprime as msp

name =  sys.argv[1]
hdf5 =  sys.argv[2]

file = open(name, 'r')
list = csv.DictReader(file, delimiter=" ")

data = msp.load(hdf5)

print("MarkerID SampleID0 Chr0 SampleID1 Chr1 Shared TrueLHS TrueRHS")

runs = []
count = 0
for line in list:
	line['mid'] = int(line['MarkerID'])
	line['id0'] = int(line['SampleID0']) * 2
	line['id1'] = int(line['SampleID1']) * 2
	line['shr'] = int(line['Shared'])

	if line['Chr0'] == 'P':
		line['id0'] = line['id0'] + 1
	if line['Chr1'] == 'P':
		line['id1'] = line['id1'] + 1
	
	line['flag'] = False
	line['ends'] = False
	line['last'] = -1
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
					run['cx'] = mut.index
	
		for run in runs:
			print("%d %s %s %s %s %d %d %d" % (run['mid'], run['SampleID0'], run['Chr0'], run['SampleID1'], run['Chr1'], run['shr'], run['xlhs'], run['xrhs']))
	
		runs = []
		count = 0

