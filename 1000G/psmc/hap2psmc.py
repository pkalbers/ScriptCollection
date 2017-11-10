import sys

dataf =  sys.argv[1]
tempf =  sys.argv[2]

lines = [line.rstrip('\n') for line in open(dataf)]

positions = []
alleles = []

for p in lines[0].split(','):
	positions.append(int(p))
	alleles.append([])

h0 = lines[1]
h1 = lines[2]

for i, letter in enumerate(h0):
	alleles[i].append(letter)

for i, letter in enumerate(h1):
	alleles[i].append(letter)

out = open(tempf, 'w')

lastPos = 0
for i in range(len(positions)):
	realPos = positions[i]
	if realPos > lastPos:
		out.write("20\t" + str(realPos) + "\t" + str(realPos - lastPos) + "\t" + "".join(alleles[i]) + "\n")
	lastPos = realPos

out.close()

