import msprime

data = msprime.load("OutOfAfricaHapMap20-v3.hdf5")

print("Mutations: ", data.get_num_mutations())
print("Samples:   ", data.get_sample_size())


file = open("OutOfAfricaHapMap20.times.txt",'w')

file.write("position node leaves parent.node node.time parent.time branch.length\n")

for tree in data.trees():
	for mut in tree.mutations():
		print(mut.position, mut.node, tree.get_num_leaves(mut.node), tree.get_parent(mut.node), tree.get_time(mut.node), tree.get_time(tree.get_parent(mut.node)), tree.get_time(tree.get_parent(mut.node)) - tree.get_time(mut.node), file=file)


