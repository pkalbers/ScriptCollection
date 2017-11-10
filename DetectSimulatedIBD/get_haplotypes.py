import msprime
import sys, getopt
from array import *
import os

num_args = 2
arg_i = ""     # input filename
arg_t = False  # transpose haplotypes
arg_b = 10000  # buffer size

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:t:b:")

for opt, arg in opts:
    if opt == '-i':
        arg_i = arg
        num_args -= 1
    elif opt == '-t':
        if arg in ["TRUE", "True", "true", "T", "t", "1"]:
            arg_t = True
        elif arg in ["FALSE", "False", "false", "F", "f", "0"]:
            arg_t = False
        else:
            print("Cannot interpret: " + arg)
            sys.exit(2)
        num_args -= 1
    elif opt == '-b':
        arg_b = int(arg)
    else:
        print("Unknown command line arguments!")
        sys.exit(2)

if num_args != 0:
    print("Missing command line arguments!")
    sys.exit(2)

arg_o = os.path.basename(arg_i)

if arg_t:
    arg_o += ".VarByRow.hap"
else:
    arg_o += ".VarByCol.hap"

arg_i += ".hdf5"


print("")
print("PRINT HAPLOTYPES")
print("")
print("Input tree sequence file: " + arg_i)
print("Output haplotypes file:   " + arg_o)
print("Transpose haplotypes:     " + str(arg_t))
if arg_t:
    print("(Buffer size: " + str(arg_b) + ")")
print("")

print("Loading data ...")
data = msprime.load(arg_i)
print("# individuals: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("")

hapfile = open(arg_o, "w")

if arg_t:
    print("Writing haplotypes (variants by line) ...")
    nsam = data.get_sample_size()
    nmut = data.get_num_mutations()
    for i in range(0, nmut, arg_b):
        print(str(i) + " of " + str(nmut))
        # allocate buffer
        buff = []
        for j in range(arg_b):
            buff.append(array("c", ['.'] * (nsam + 1)))
        # walkabout samples
        k = 0
        for hap in data.haplotypes():
            for j in range(arg_b):
                ij = i + j
                if ij == nmut:
                    break
                buff[j][k] = hap[ij]
            k += 1
        # write buffer to file
        for j in range(arg_b):
            if buff[j][0] != ".":
                buff[j][nsam] = '\n'
                hapfile.write(buff[j].tostring())
else:
    print("Writing haplotypes (individuals by line) ...")
    for hap in data.haplotypes():
        hap += "\n"
        hapfile.write(hap)

hapfile.close()

print("Complete!")
