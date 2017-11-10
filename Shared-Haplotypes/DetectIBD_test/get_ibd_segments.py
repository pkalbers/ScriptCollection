import msprime
import sys, getopt

num_args = 3
arg_i = ""    # input filename
arg_p = ""    # pairs filename
arg_o = ""    # output filename

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:o:p:")

for opt, arg in opts:
    if opt == '-i':
        arg_i = arg
        num_args -= 1
    elif opt == '-o':
        arg_o = arg
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


print("")
print("Input tree sequence file: " + arg_i)
print("Sharing pairs table:      " + arg_p)
print("Output segments file:     " + arg_o)
print("")


print("Loading data ...")
data = msprime.load(arg_i)
print("# individuals: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("")


print("")
print("Reading sharing pairs table ...")

pairs = list()
with open(arg_p) as pfile:
    for pline in pfile:
        idx, pos, p0, p1 = pline.split()
        pairs.append([int(idx), float(pos), int(p0), int(p1), int(), float(), float(), float(), float(), float()])

print("Done")
print("")


print("")
print("Finding the MRCA for each pair ...")

c = 0
p = 0
n = len(pairs)
f = False
for tree in data.trees():
    if (f):
        break
    for pos, node in tree.mutations():
        pos = round(float(str(pos)), 4)
        if (f):
            break
        if (pos < round(pairs[p][1], 4)):
            continue
        while (pos > round(pairs[p][1], 4)):
            print("Mismatch between positions " + str(pos) + " and " + str(round(pairs[p][1], 4)))
            p += 1
            if (p == n):
                f = True
                break
        if (f):
            break
        while (pos == round(pairs[p][1], 4)):
            pairs[p][4] = tree.get_mrca(pairs[p][2], pairs[p][3])   #  MRCA
            pairs[p][5] = tree.get_tmrca(pairs[p][2], pairs[p][3])  # tMRCA
            c += 1
            if (c % 10000 == 0):
                print(str(c) + " of " + str(n))
            p += 1
            if (p == n):
                f = True
                break

print("Done")
print("")


print("")
print("IBD detection for MRCA for each pair ...")

n = len(pairs)
for p in range(n):
    if (p % 10000 == 0):
        print(str(p) + " of " + str(n))
    ppos = round(pairs[p][1], 4)
    pair = [ pairs[p][2], pairs[p][3] ]
    mrca = pairs[p][4]
    tmp = [ float(-1), float(-1) ]
    lhs = True
    rhs = False
    ibd = False
    for tree in data.trees():
        if (rhs == True):
            break
        is_mrca = (mrca == tree.get_mrca(pair[0], pair[1]))
        for tpos, node in tree.mutations():
            if (lhs == True):
                if (ibd == False and is_mrca == True):
                    tmp[0] = tpos
                    ibd = True
                if (ibd == True and is_mrca == False):
                    tmp[0] = float(-1)
                    ibd = False
                if (ppos == round(float(str(tpos)), 4)):
                    lhs = False
                    tmp[1] = tpos
            else:
                if (ibd == False):
                    print("Data mismatch!!!")
                    rhs = True
                    break
                if (is_mrca == True):
                    tmp[1] = tpos
                else:
                    rhs = True
                    break
    if (tmp[0] == float(-1) or tmp[1] == float(-1)):
        print("IBD not detected: " + str(pairs[p]))
    pairs[p][6] = tmp[0]
    pairs[p][7] = tmp[1]

print("Done")
print("")

print("")
print("Saving ...")
ibdfile = open(arg_o, "w")
ibdfile.write("index position h0 h1 mrca tmrca mrca.lhs mrca.rhs tmrca.lhs tmrca.rhs\n")
for p in range(len(pairs)):
	ibdfile.write(" ".join(map(str, pairs[p])) + "\n")
ibdfile.close()
print("")

print("")
print("IBD detection for tMRCA for each pair ...")

n = len(pairs)
for p in range(n):
    if (p % 10000 == 0):
        print(str(p) + " of " + str(n))
    ppos = round(pairs[p][1], 4)
    pair = [ pairs[p][2], pairs[p][3] ]
    tmrca = pairs[p][5]
    tmp = [ float(-1), float(-1) ]
    lhs = True
    rhs = False
    ibd = False
    for tree in data.trees():
        if (rhs == True):
            break
        is_tmrca = (tmrca == tree.get_tmrca(pair[0], pair[1]))
        for tpos, node in tree.mutations():
            if (lhs == True):
                if (ibd == False and is_tmrca == True):
                    tmp[0] = tpos
                    ibd = True
                if (ibd == True and is_tmrca == False):
                    tmp[0] = float(-1)
                    ibd = False
                if (ppos == round(float(str(tpos)), 4)):
                    lhs = False
                    tmp[1] = tpos
            else:
                if (ibd == False):
                    print("Data mismatch!!!")
                    rhs = True
                    break
                if (is_tmrca == True):
                    tmp[1] = tpos
                else:
                    rhs = True
                    break
    if (tmp[0] == float(-1) or tmp[1] == float(-1)):
        print("IBD not detected: " + str(pairs[p]))
    pairs[p][8] = tmp[0]
    pairs[p][9] = tmp[1]

print("Done")
print("")


print("")
print("Saving ...")
ibdfile = open(arg_o, "w")
ibdfile.write("index position h0 h1 mrca tmrca mrca.lhs mrca.rhs tmrca.lhs tmrca.rhs\n")
for p in range(len(pairs)):
	ibdfile.write(" ".join(map(str, pairs[p])) + "\n")
ibdfile.close()
print("")


print("Complete!")
