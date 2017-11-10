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

arg_o = os.path.basename(arg_p)
arg_o = "age." + arg_o


print("")
print("DETECT IBD")
print("")
print("Input tree sequence file: " + arg_i)
print("Sharing pairs table:      " + arg_p)
print("Output file:              " + arg_o)
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
        pair.append([int(idx), float(pos), int(fk), int(h0), int(h1), float(-1)])

pair_size = len(pair)

print("Done")
print("")


print("")
print("Finding tMRCA for each pair ...")

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
            pair[p][5] = tree.get_tmrca(pair[p][3], pair[p][4])  # tMRCA
            c += 1
            if (c % 1000 == 0):
                print(str(c) + " of " + str(pair_size))
            p += 1
            if (p == pair_size):
                f = True
                break

print("Done")
print("")


print("")
print("Saving ...")
agefile = open(arg_o, "w")
agefile.write("index position fk h0 h1 age\n")
for p in range(pair_size):
	agefile.write(" ".join(map(str, pair[p])) + "\n")
agefile.close()
print("")

