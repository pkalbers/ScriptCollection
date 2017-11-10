import msprime
import sys, getopt, random, math

num_args = 2
arg_o = ""  # output filename prefix
arg_r = ""  # recombination map, filename

# Read command line args
opts, args = getopt.getopt(sys.argv[1:], "o:r:")

for opt, arg in opts:
    if opt == '-o':
        arg_o = arg
        num_args -= 1
    elif opt == '-r':
        arg_r = arg
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
print("Recombination Map:         " + arg_r)
print("Output file:               " + arg_o)
print("")


print("Reading recombination map ... ")
recomb_map = msprime.RecombinationMap.read_hapmap(arg_r)
print("OK")
print("")


# mutation rate according to Gutenkunst et al 2009
MUTRATE = 2.35e-8

# First we set out the maximum likelihood values of the various parameters
# given in Table 1.
N_A = 7300
N_B = 2100
N_AF = 12300
N_EU0 = 1000
N_AS0 = 510
# Times are provided in years, so we convert into generations.
generation_time = 25
T_AF = 220e3 / generation_time
T_B = 140e3 / generation_time
T_EU_AS = 21.2e3 / generation_time
# We need to work out the starting population sizes based on the growth
# rates provided for these two populations
r_EU = 0.004
r_AS = 0.0055
N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
# Migration rates during the various epochs.
m_AF_B = 25e-5
m_AF_EU = 3e-5
m_AF_AS = 1.9e-5
m_EU_AS = 9.6e-5
# Population IDs correspond to their indexes in the popupulation
# configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
# initially.
population_configurations = [
	msprime.PopulationConfiguration(sample_size=0, initial_size=N_AF),
	msprime.PopulationConfiguration(sample_size=5000, initial_size=N_EU, growth_rate=r_EU),
	msprime.PopulationConfiguration(sample_size=0, initial_size=N_AS, growth_rate=r_AS)
]
migration_matrix = [
	[      0, m_AF_EU, m_AF_AS],
	[m_AF_EU,       0, m_EU_AS],
	[m_AF_AS, m_EU_AS,       0],
]
demographic_events = [
	# CEU and CHB merge into B with rate changes at T_EU_AS
	msprime.MassMigration(
		time=T_EU_AS, source=2, destination=1, proportion=1.0),
	msprime.MigrationRateChange(time=T_EU_AS, rate=0),
	msprime.MigrationRateChange(
		time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
	msprime.MigrationRateChange(
		time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
	msprime.PopulationParametersChange(
		time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
	# Population B merges into YRI at T_B
	msprime.MassMigration(
		time=T_B, source=1, destination=0, proportion=1.0),
	# Size changes to N_A at T_AF
	msprime.PopulationParametersChange(
		time=T_AF, initial_size=N_A, population_id=0)
]

print("Simulating tree sequence ... ")

random.seed(random.randrange(100))
seed = int(random.random() * random.randrange(10, 1000000))
#data = msprime.simulate(sample_size=arg_s, Ne=arg_n, mutation_rate=arg_m, recombination_map=recomb_map, random_seed=seed)

data = msprime.simulate(Ne=N_A, mutation_rate=MUTRATE,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        recombination_map=recomb_map,
        random_seed=seed)

print("Saving dump file ...")

data.dump(arg_o)

print("Complete!")
