import msprime
import sys

hdf5 =  sys.argv[1]
temp =  sys.argv[2]

i0 =  int(sys.argv[3])
i1 =  int(sys.argv[4])

data = msprime.load(hdf5)

positions = []
alleles = []

for mut in data.mutations():
	positions.append(int(mut.position))
	alleles.append([])

x0 = False
x1 = False

h0 = ""
h1 = ""

i = 0
for hap in data.haplotypes():
	if i == i0:
		h0 = hap
		x0 = True
	if i == i1:
		h1 = hap
		x1 = True
	if x0 and x1:
		break
	i += 1

if (x0 == False) or (x1 == False):
	quit()

out = open(temp, 'w')

for i, letter in enumerate(h0):
	alleles[i].append(letter)
for i, letter in enumerate(h1):
	alleles[i].append(letter)

lastPos = 0
for i in range(len(positions)):
	realPos = positions[i]
	if realPos > lastPos:
		out.write("1\t" + str(realPos) + "\t" + str(realPos - lastPos) + "\t" + "".join(alleles[i]) + "\n")
	lastPos = realPos

out.close()

