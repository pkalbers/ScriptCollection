import msprime
import sys, getopt

num_args = 2
arg_i = ""    # input filename
arg_o = ""    # output filename

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:o:")

for opt, arg in opts:
    if opt == '-i':
        arg_i = arg
        num_args -= 1
    elif opt == '-o':
        arg_o = arg
        num_args -= 1
    else:
        print("Unknown command line arguments!")
        sys.exit(2)

if num_args != 0:
    print("Missing command line arguments!")
    sys.exit(2)

arg_i += ".hdf5"
arg_o += ".pos"

print("")
print("Input tree sequence file: " + arg_i)
print("Output positions file:    " + arg_o)
print("")

print("Loading data ...")
data = msprime.load(arg_i)
print("# individuals: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("")

posfile = open(arg_o, "w")

i = 0
for tree in data.trees():
    for position, node in tree.mutations():
        i += 1
        posfile.write(str(position) + "\n")

print(str(i) + " positions written to " + arg_o)
