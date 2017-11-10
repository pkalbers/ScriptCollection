import msprime
import sys, getopt
import os

num_args = 3
arg_i = ""    # input filename
arg_p = ""    # pairs filename
arg_l = 0.0   # minimum length

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:p:l:")

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
print("# chromosomes: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("# trees:       " + str(data.get_num_trees()))
print("")


n_tree = data.get_num_trees()


print("Reading sharing pairs table ...")

G = list()
H = list()
n_pair = 0
with open(arg_p) as pfile:
    for pline in pfile:
        g0, g1, h0, h1 = pline.split()
        G.append( ( int(g0), int(g1) ) )
        H.append( ( int(h0), int(h1) ) )
        n_pair += 1

r_pair = range(n_pair)

print("Number of pairs: " + str(n_pair))
print("")



print("Global IBD detection for each pair ...")

segment = list()

mrca = list()
for i in r_pair:
	mrca.append([ float(-1) , float(-1), int(-1) ])  # beg, end, mrca

count_tree = 0
count_all = 0
count_sub = 0
min_length = float(arg_l)
for tree in data.trees():
	intv = tree.get_interval()
	if (count_tree == 0):
		for i in r_pair:
			mrca[i][0] = intv[0]
			mrca[i][1] = intv[1]
			mrca[i][2] = tree.get_mrca(H[i][0], H[i][1])
		count_tree += 1
		continue
	count_tree += 1
	if (count_tree % 1000 == 0):
		print(" " + str(count_tree) + " of " + str(n_tree) + " trees ({0:.1f}%), ".format(float(count_tree) / float(n_tree) * 100) + str(count_sub) + " segments collected, " + str(count_all) + " segments found")
	for i in r_pair:
		h = H[i]
		target = mrca[i]
		this_mrca = tree.get_mrca(h[0], h[1])
		if (count_tree == n_tree):
			target[1] = intv[1]
		elif (target[2] == this_mrca):
			mrca[i][1] = intv[1]
			continue
		if ((target[1] - target[0]) >= min_length):
			count_sub += 1
			segment.append( {
			'G' : G[i],
			'H' : H[i],
			'interval' : ( target[0] , target[1] ) } )
		count_all += 1
		mrca[i] = [ intv[0] , intv[1] , this_mrca ]

n_segment = len(segment)
r_segment = range(n_segment)

print("Number of IBD segments: " + str(n_segment))
print("")



print("Positioning segments ...")

position = list()
for tree in data.trees():
	for pos, node in tree.mutations():
		position.append(pos)

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
	if (breakpoint[brk] > position[pos]):
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

count = 0
fails = 0
for i in r_segment:
	count += 1
	if (count % 100000 == 0):
		print(" " + str(count) + " of " + str(n_segment) + " segments ({0:.1f}%)".format(float(count) / float(n_segment) * 100))
	tmp = segment[i]['interval']
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

print("Number of positioned IBD segments: " + str(count))
print("")



print("Writing to file ...")

out_stream = open(out_file, "w")
out_stream.write("g0 g1 h0 h1 beg end lhs rhs\n")

for i in r_segment:
	seg = segment[i]
	loc = location[i]
	if (loc[0] == int(-1)) or (loc[1] == int(-1)):
		continue
	out_stream.write("%d %d %d %d %.4f %.4f %d %d\n" % (
	seg['G'][0], seg['G'][1],
	seg['H'][0], seg['H'][1],
	seg['interval'][0], seg['interval'][1],
	loc[0], loc[1] ) )

out_stream.close()

print("Done")
print("")

print("Complete!")
