import msprime
import sys, getopt
import os

num_args = 3
arg_i = ""    # input filename
arg_p = ""    # pairs filename
arg_l = 0.0     # minimum length

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:p:")

for opt, arg in opts:
    if opt == '-i':
        arg_i = arg
        num_args -= 1
    elif opt == '-p':
        arg_p = arg
        num_args -= 1
    elif opt == '-l':
        arg_l = int(arg)
        num_args -= 1
    else:
        print("Unknown command line arguments!")
        sys.exit(2)

if num_args != 0:
    print("Missing command line arguments!")
    sys.exit(2)


out_file = os.path.basename(arg_p)
out_file = "global_ibd." + out_file


print("")
print("DETECT GLOBAL IBD")
print("")
print("Input history file:  " + arg_i)
print("Sharing pairs table: " + arg_p)
print("Minimum length:      " + str(arg_l))
print("Output file:         " + out_file)
print("")


print("Loading data ...")
data = msprime.load(arg_i)
print("# individuals: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("# trees:       " + str(data.get_num_trees()))
print("")


position = list()
for tree in data.trees():
	for pos, node in tree.mutations():
		position.append(pos)


print("Reading sharing pairs table ...")

pair = list()
with open(arg_p) as pfile:
    for pline in pfile:
        g0, g1, h0, h1 = pline.split()
        pair.append( {
        'g0' : int(g0), 'g1' : int(g1),
        'h0' : int(h0), 'h1' : int(h1),
        'mrca' : int(-1),
        'beg' : float(-1),
        'end' : float(-1) } )

print("Number of pairs: " + str(len(pair)))
print("")


print("Global IBD detection for each pair ...")

ibd = list()

j = 0
count_all = 0
count_sub = 0
min = float(arg_l)
max = data.get_num_trees()
for tree in data.trees():
	if (j % 1000 == 0):
		print(" " + str(j) + " of " + str(data.get_num_trees()) + " trees, " + str(count_sub) + " segments collected, " + str(count_all) + " segments found")
	intv = tree.get_interval()
	beg = intv[0]
	end = intv[1]
	if (j == 0):
		for i in range(len(pair)):
			pair[i]['mrca'] = tree.get_mrca(pair[i]['h0'], pair[i]['h1'])
			pair[i]['beg'] = beg
			pair[i]['end'] = end
	j += 1
	for i in range(len(pair)):
		store = False
		tmp = pair[i]
		last_mrca = tmp['mrca']
		this_mrca = tree.get_mrca(tmp['h0'], tmp['h1'])
		if (last_mrca == this_mrca):
			pair[i]['end'] = end
		else:
			store = True
		if (j == max):
			pair[i]['end'] = end
			store = True
		if (store == True):
			if ((pair[i]['end'] - tmp['beg']) >= min):
				count_sub += 1
				ibd.append( {
				'g0' : tmp['g0'], 'g1' : tmp['g1'],
				'h0' : tmp['h0'], 'h1' : tmp['h1'],
				'lhs' : int(-1),    'rhs' : int(-1),
				'beg' : tmp['beg'], 'end' : pair[i]['end'] } )
			count_all += 1
			pair[i]['mrca'] = this_mrca
			pair[i]['beg'] = beg
			pair[i]['end'] = end

print("Done")
print("")



print("Positioning segments ...")

for i in range(len(ibd)):
	if (i % 10000 == 0):
		print(" " + str(i) + " of " + str(len(ibd)) + " segments")
	beg = ibd[i]['beg']
	end = ibd[i]['end']
	switcher = False
	last = int(-1)
	for j in range(len(position)):
		if (switcher == False):
			if (beg <= position[j]):
				ibd[i]['lhs'] = last = j
				switcher = True
		else:
			if (end <= position[j]):
				break
			last = j
	ibd[i]['rhs'] = last

error = False
for i in range(len(ibd)):
    if (ibd[i]['lhs'] == int(-1)) or (ibd[i]['rhs'] == int(-1)):
        print("Segment not positioned at line " + str(i + 1))
        error = True

if (error == True):
    sys.exit("ABORTED")

print("Done")
print("")



print("Writing to file ...")

out_stream = open(out_file, "w")
out_stream.write("g0 g1 h0 h1 lhs rhs beg end\n")

for i in range(len(ibd)):
	out_stream.write("%d %d %d %d %d %d %.4f %.4f\n" % (
	ibd[i]['g0'], ibd[i]['g1'],
	ibd[i]['h0'], ibd[i]['h1'],
	ibd[i]['lhs'], ibd[i]['rhs'],
	ibd[i]['beg'], ibd[i]['end']))

out_stream.close()

print("Done")
print("")

print("Complete!")
