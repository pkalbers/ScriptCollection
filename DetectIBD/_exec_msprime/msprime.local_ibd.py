import msprime
import numpy as np
import sys, getopt
import os

num_args = 2
arg_i = ""    # input filename
arg_p = ""    # pairs filename

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:p:")

for opt, arg in opts:
    if opt == '-i':
        arg_i = arg
        num_args -= 1
    elif opt == '-p':
        arg_p = arg
        num_args -= 1
    else:
        print("Unknown command line arguments!")
        sys.exit(2)

if num_args != 0:
    print("Missing command line arguments!")
    sys.exit(2)


out_file = os.path.basename(arg_p)
out_file = "local_ibd." + out_file


print("")
print("DETECT LOCAL IBD")
print("")
print("Input history file:  " + arg_i)
print("Sharing pairs table: " + arg_p)
print("Output file:         " + out_file)
print("")


print("Loading data ...")
data = msprime.load(arg_i)
print("# chromosomes: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("# trees:       " + str(data.get_num_trees()))
print("")


n_tree = data.get_num_trees()


print("Reading sharing pairs table ...")

pair = list()
G = list()
H = list()
with open(arg_p) as pfile:
    for pline in pfile:
        idx, pos, fk, g0, g1, h0, h1 = pline.split()
        pair.append( { 'idx' : int(idx), 'pos' : pos, 'fk' : int(fk) } )
        G.append( ( int(g0), int(g1) ) )
        H.append( ( int(h0) - 1, int(h1) - 1 ) )  ### conversion from R where index begins at 1

n_pair = len(pair)
r_pair = range(n_pair)

print("Number of pairs: " + str(n_pair))
print("")



print("Finding focal MRCA for each pair ...")

position = list()
for tree in data.trees():
	for pos, node in tree.mutations():
		position.append(pos)

focal = list()
focal_mrca = list()
focal_time = list()
for i in r_pair:
	focal.append(position[ pair[i]['idx'] - 1 ])
	focal_mrca.append( int(-1) )
	focal_time.append( float(-1) )

focal_array = np.array(focal)

count_tree = 0
count_pair = 0
for tree in data.trees():
	if (count_pair == n_pair):
		print("All found")
		break
	intv = tree.get_interval()
	count_tree += 1
	if (count_tree % 10000 == 0):
		print(" " + str(count_tree) + " of " + str(n_tree) + " trees ({0:.1f}%), ".format(float(count_tree) / float(n_tree) * 100) + str(count_pair) + " MRCA found")
	i = (intv[0] <= focal_array) & (focal_array < intv[1])
	i = np.flatnonzero(i)
	if (i.size == 0):
		continue
	for j in np.nditer(i):
		h = H[j]
		this_mrca = tree.get_mrca(h[0], h[1])
		focal_mrca[j] = this_mrca
		focal_time[j] = tree.get_time(this_mrca)  # time to MRCA
		count_pair += 1

print("Focal MRCA detected: " + str(count_pair) + " of " + str(n_pair) + " pairs")

error = False
for i in r_pair:
    if (focal_mrca[i] == int(-1)) or (focal_time[i] == float(-1)):
        print("No MRCA found at" + str(i) + " for " + str(pair[i]['idx']))
        error = True

if (error == True) or not (count_pair == n_pair):
    sys.exit("ABORTED")

print("")



print("Local IBD detection for each pair ...")

segment = list()

ibd_state = list()
ibd_found = list()
for i in r_pair:
	segment.append([ float(-1), float(-1) ])
	ibd_state.append(False)
	ibd_found.append(False)

par = list()
for i in r_pair:
	par.append( [ H[i][0], H[i][1], focal[i], focal_mrca[i] ] )

count_tree = 0
count_pair = 0
count_found = 0
for tree in data.trees():
	if (count_found == n_pair):
		print("All found")
		break
	intv = tree.get_interval()
	count_tree += 1
	if (count_tree % 1000 == 0):
		print(" " + str(count_tree) + " of " + str(n_tree) + " trees ({0:.1f}%), ".format(float(count_tree) / float(n_tree) * 100) + str(count_pair) + " segments found")
	for i in r_pair:
		if (ibd_found[i] == True):
			continue
		h = H[i]
		find_pos  = focal[i]
		find_mrca = focal_mrca[i]
		this_mrca = tree.get_mrca(h[0], h[1])
		is_lhs = (intv[0] <= find_pos)
		is_rhs = (intv[1] >  find_pos)
		if (this_mrca == find_mrca):
			if (is_lhs == True) and (is_rhs == True):
				count_pair += 1
			if (ibd_state[i] == False):
				ibd_state[i] = True
				segment[i][0] = intv[0]
				segment[i][1] = intv[1]
			else:
				segment[i][1] = intv[1]
		else:
			if (is_lhs == True) and (is_rhs == True):
				sys.exit("MRCA mismatch!")
			if (ibd_state[i] == True):
				ibd_state[i] = False
				if is_lhs:
					segment[i][0] = float(-1)
					segment[i][1] = float(-1)
				if is_rhs:
					ibd_found[i] = True
					count_found += 1

print("IBD segments detected: " + str(count_pair) + " of " + str(n_pair) + " pairs")

error = False

for i in r_pair:
    if (segment[i][0] == float(-1)) or (segment[i][1] == float(-1)):
        print("No segment found at" + str(i) + " for " + str(pair[i]['idx']))
        error = True
    if (segment[i][0] > segment[i][1]):
        print("Reversed segment interval at" + str(i) + " for " + str(pair[i]['idx']))
        error = True
    if (segment[i][0] > focal[i]) or (segment[i][1] <= focal[i]):
        print("Focal site out of segment bounds at" + str(i) + " for " + str(pair[i]['idx']))
        error = True

if (error == True) or not (count_pair == n_pair):
    sys.exit("ABORTED")

print("")



print("Positioning segments ...")

breakpoint = list()
for brk in data.breakpoints():
	breakpoint.append(brk)

lhs_break = dict()
rhs_break = dict()

brk = 0
pos = 0
max_pos = len(position)
max_brk = len(breakpoint)
last = pos
while (pos < max_pos) and (brk < max_brk):
	if (position[pos] < breakpoint[brk]):
		last = pos
		pos += 1
	else:
		key = "%.8f" % breakpoint[brk]
		lhs_break[key] = pos
		rhs_break[key] = last
		brk += 1

while (brk < max_brk):
	key = "%.8f" % breakpoint[brk]
	rhs_break[key] = last
	brk += 1

location = list()

fails = 0
for i in r_pair:
	tmp = segment[i]
	beg = "%.8f" % tmp[0]
	end = "%.8f" % tmp[1]
	lhs = int(-1)
	rhs = int(-1)
	try:
		lhs = lhs_break[beg]
		rhs = rhs_break[end]
	except:
		fails += 1
	location.append( ( lhs , rhs ) )

if (fails > 0):
	print("# failed segments: " + str(fails))

error = False
for i in r_pair:
    if (location[i][0] == int(-1)) or (location[i][1] == int(-1)):
        print("Segment not positioned at " + str(i) + " for " + str(pair[i]['idx']))
        error = True
    if (location[i][0] > pair[i]['idx'] - 1) or (location[i][1] < pair[i]['idx'] - 1):
        print("Segment wrongly positioned at " + str(i) + " for " + str(pair[i]['idx']))
        error = True

if (error == True):
    sys.exit("ABORTED")

print("Done")
print("")



print("Writing to file ...")

out_stream = open(out_file, "w")
out_stream.write("index position fk g0 g1 h0 h1 time lhs.position rhs.position lhs.index rhs.index\n")

for i in r_pair:
	seg = segment[i]
	loc = location[i]
	out_stream.write("%d %s %d %d %d %d %d %.16f %.16f %.16f %d %d\n" % (
	pair[i]['idx'],
	pair[i]['pos'],
	pair[i]['fk'],
	G[i][0], G[i][1],
	H[i][0] + 1, H[i][1] + 1,  ### back-conversion to R
	focal_time[i],
	seg[0], seg[1],
	loc[0], loc[1] ) )

out_stream.close()

print("Done")
print("")

print("Complete!")
