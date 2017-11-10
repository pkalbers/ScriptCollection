import msprime

data = msprime.simulate(Ne = 10000,
						mutation_rate = 1e-08,
						sample_size = 2000,
						recombination_rate = 1e-08,
						length = 100000000)

print("Mutations: ", data.get_num_mutations())
print("Samples:   ", data.get_sample_size())
print("Breaks:    ", len(list(data.breakpoints())))


data.dump("vanilla.hdf5")

with open("vanilla.mutations.txt", "w") as mutations_file:
	data.write_mutations(mutations_file)

with open("vanilla.records.txt", "w") as records_file:
	data.write_records(records_file)

with open("vanilla.breakpoints.txt", "w") as brk_file:
	for brk in list(data.breakpoints()):
		brk_file.write("{}\n".format(brk))

with open("vanilla.vcf", "w") as vcf_file:
	data.write_vcf(vcf_file, 2)



