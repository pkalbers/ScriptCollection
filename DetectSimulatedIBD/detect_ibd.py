import msprime
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


arg_i += ".hdf5"

file_mrca = os.path.basename(arg_p)
file_mrca = "ibd_mrca." + file_mrca

file_tmrca = os.path.basename(arg_p)
file_tmrca = "ibd_tmrca." + file_tmrca


print("")
print("DETECT IBD")
print("")
print("Input tree sequence file: " + arg_i)
print("Sharing pairs table:      " + arg_p)
print("Output files:             " + file_mrca)
print("                          " + file_tmrca)
print("")


print("Loading data ...")
data = msprime.load(arg_i)
print("# individuals: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("")


print("")
print("Reading sharing pairs table ...")

pair = list()
with open(arg_p) as pfile:
    for pline in pfile:
        idx, pos, fk, h0, h1 = pline.split()
        pair.append([int(idx), float(pos), int(fk), int(h0), int(h1)])

pair_size = len(pair)

print("Done")
print("")


print("")
print("Finding the MRCA + tMRCA for each pair ...")

mrca = list()
for p in range(pair_size):
    mrca.append(int(-1))

tmrca = list()
for p in range(pair_size):
    tmrca.append(float(-1))

c = 0
p = 0
f = False
for tree in data.trees():
    if (f):
        break
    for pos, node in tree.mutations():
        pos = "{:.8f}".format(pos)
        pos = round(float(pos), 4)
        if (f):
            break
        if (pos < round(pair[p][1], 4)):
            continue
        while (pos > round(pair[p][1], 4)):
            print("Mismatch between positions " + str(pos) + " and " + str(round(pair[p][1], 4)))
            p += 1
            if (p == pair_size):
                f = True
                break
        if (f):
            break
        while (pos == round(pair[p][1], 4)):
            mrca[p]  = tree.get_mrca(pair[p][3],  pair[p][4])  #  MRCA
            tmrca[p] = tree.get_tmrca(pair[p][3], pair[p][4])  # tMRCA
            c += 1
            if (c % 1000 == 0):
                print(str(c) + " of " + str(pair_size))
            p += 1
            if (p == pair_size):
                f = True
                break

print("Done")
print("")


error = False
for p in range(pair_size):
    if (mrca[p] == int(-1)):
        print("No MRCA found for " + str(pair[p]))
        error = True
if (error == True):
    sys.exit("ABORTED")

error = False
for p in range(pair_size):
    if (tmrca[p] == float(-1)):
        print("No tMRCA found for " + str(pair[p]))
        error = True
if (error == True):
    sys.exit("ABORTED")


print("")
print("IBD detection for MRCA for each pair ...")

stream_mrca = open(file_mrca, "w")
stream_mrca.write("index position fk h0 h1 lhs rhs\n")

for p in range(pair_size):
    if (p % 1000 == 0):
        print(str(p) + " of " + str(pair_size))
    tmp = pair[p]
    pair_pos = round(tmp[1], 4)
    pair_h0 = tmp[3]
    pair_h1 = tmp[4]
    pair_lhs = float(-1)
    pair_rhs = float(-1)
    pair_mrca = mrca[p]
    flag_lhs = True
    flag_rhs = False
    flag_ibd = False
    for tree in data.trees():
        if (flag_rhs == True):
            break
        is_mrca = (pair_mrca == tree.get_mrca(pair_h0, pair_h1))
        for tree_pos, node in tree.mutations():
            if (flag_lhs == True):
                if (flag_ibd == False and is_mrca == True):
                    pair_lhs = tree_pos
                    flag_ibd = True
                if (flag_ibd == True and is_mrca == False):
                    pair_lhs = float(-1)
                    flag_ibd = False
                if (pair_pos == round(float("{:.8f}".format(tree_pos)), 4)):
                    flag_lhs = False
                    pair_rhs = tree_pos
            else:
                if (flag_ibd == False):
                    flag_rhs = True
                    break
                if (is_mrca == True):
                    pair_rhs = tree_pos
                else:
                    flag_rhs = True
                    break
    if (pair_lhs == float(-1) or pair_rhs == float(-1)):
        print("IBD not detected: " + str(tmp))
    else:
        stream_mrca.write(" ".join(map(str, tmp)) + " " + str(pair_lhs) + " " + str(pair_rhs) + "\n")

stream_mrca.close()
print("Done")
print("")


sys.exit("Skipping tMRCA detection")


print("")
print("IBD detection for tMRCA for each pair ...")

stream_tmrca = open(file_tmrca, "w")
stream_tmrca.write("index position fk h0 h1 lhs rhs\n")

for p in range(pair_size):
    if (p % 1000 == 0):
        print(str(p) + " of " + str(pair_size))
    tmp = pair[p]
    pair_pos = round(tmp[1], 4)
    pair_h0 = tmp[3]
    pair_h1 = tmp[4]
    pair_lhs = float(-1)
    pair_rhs = float(-1)
    pair_tmrca = tmrca[p]
    flag_lhs = True
    flag_rhs = False
    flag_ibd = False
    for tree in data.trees():
        if (flag_rhs == True):
            break
        is_tmrca = (pair_tmrca == tree.get_tmrca(pair_h0, pair_h1))
        for tree_pos, node in tree.mutations():
            if (flag_lhs == True):
                if (flag_ibd == False and is_tmrca == True):
                    pair_lhs = tree_pos
                    flag_ibd = True
                if (flag_ibd == True and is_tmrca == False):
                    pair_lhs = float(-1)
                    flag_ibd = False
                if (pair_pos == round(float("{:.8f}".format(tree_pos)), 4)):
                    flag_lhs = False
                    pair_rhs = tree_pos
            else:
                if (flag_ibd == False):
                    flag_rhs = True
                    break
                if (is_tmrca == True):
                    pair_rhs = tree_pos
                else:
                    flag_rhs = True
                    break
    if (pair_lhs == float(-1) or pair_rhs == float(-1)):
        print("IBD not detected: " + str(tmp))
    else:
        stream_tmrca.write(" ".join(map(str, tmp)) + " " + str(pair_lhs) + " " + str(pair_rhs) + "\n")

stream_tmrca.close()
print("Done")
print("")


print("Complete!")
