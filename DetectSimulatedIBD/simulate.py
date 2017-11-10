import msprime
import sys, getopt, random

num_args = 5
arg_o = ""  # output filename prefix
arg_s = 0   # sample size
arg_n = 0   # number of loci
arg_r = 0   # recombination rate, scaled
arg_m = 0   # mutation rate, scaled

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "o:s:n:r:m:")

for opt, arg in opts:
    if opt == '-o':
        arg_o = arg
        num_args -= 1
    elif opt == '-s':
        arg_s = int(arg)
        num_args -= 1
    elif opt == '-n':
        arg_n = int(arg)
        num_args -= 1
    elif opt == '-r':
        arg_r = float(arg)
        num_args -= 1
    elif opt == '-m':
        arg_m = float(arg)
        num_args -= 1
    else:
        print("Unknown command line arguments!")
        sys.exit(2)

if num_args != 0:
    print("Missing command line arguments!")
    sys.exit(2)

arg_o += ".hdf5"

print("")
print("SIMULATE")
print("")
print("Sample size:               " + str(arg_s))
print("Number of loci:            " + str(arg_n))
print("Scaled recombination rate: " + str(arg_r))
print("Scaled mutation rate:      " + str(arg_m))
print("Output file:               " + arg_o)
print("")

print("Simulating tree sequence ...")

random.seed(random.randrange(100))
seed = int(random.random() * random.randrange(10, 1000000))
data = msprime.simulate(arg_s, arg_n, arg_r, arg_m, random_seed=seed)

print("Saving dump file ...")

data.dump(arg_o)

print("Complete!")
