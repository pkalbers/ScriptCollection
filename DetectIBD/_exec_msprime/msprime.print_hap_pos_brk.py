import msprime
import sys, getopt
from array import *
import os

num_args = 1
arg_i = ""     # input filename

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "i:")

for opt, arg in opts:
    if opt == '-i':
        arg_i = arg
        num_args -= 1
    else:
        print("Unknown command line arguments!")
        sys.exit(2)

if num_args != 0:
    print("Missing command line arguments!")
    sys.exit(2)


prefix = os.path.basename(arg_i)

hap_file = prefix + ".hap"
pos_file = prefix + ".pos"
brk_file = prefix + ".brk"


print("")
print("Input history file: " + arg_i)
print("Haplotypes file:    " + hap_file)
print("Positions file:     " + pos_file)
print("Breakpoints file:   " + brk_file)
print("")

print("Loading data ...")
data = msprime.load(arg_i)
print("# individuals: " + str(data.get_sample_size()))
print("# variants:    " + str(data.get_num_mutations()))
print("# trees:       " + str(data.get_num_trees()))
print("")

print("> Writing haplotypes")
hap_file = open(hap_file, "w")
for hap in data.haplotypes():
    hap += "\n"
    hap_file.write(hap)
hap_file.close()

print("> Writing positions")
pos_file = open(pos_file, "w")
for tree in data.trees():
    for pos, tmp in tree.mutations():
        pos = '%.16f' % pos  # '%s' % float('%.16f' % pos)
        pos_file.write(pos + "\n")
pos_file.close()

print("> Writing breakpoints")
brk_file = open(brk_file, "w")
for brk in data.breakpoints():
    brk = '%.16f' % brk  # '%s' % float('%.16f' % brk)
    brk_file.write(brk + "\n")
brk_file.close()

print("Complete!")
