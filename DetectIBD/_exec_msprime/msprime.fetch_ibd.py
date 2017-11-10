import msprime

### assign:
# --> history
# --> h0
# --> h1
###

def MRCA(h0, h1, history):
	mrca = list()
	time = list()
	lhs_position = list()
	rhs_position = list()
	lhs_index = list()
	rhs_index = list()
	data = msprime.load(history)
	beg = 1
	end = 1
	i = 0
	for tree in data.trees():
		this_mrca = tree.get_mrca(h0 - 1, h1 - 1)
		this_time = tree.get_time(this_mrca)
		this_intv = tree.get_interval()
		mrca.append(this_mrca)
		time.append(this_time)
		lhs_position.append(this_intv[0])
		rhs_position.append(this_intv[1])
		j = 0
		for pos, node in tree.mutations():
			i += 1
			if (j == 0):
				beg = i
			end = i
			j += 1
		lhs_index.append(beg)
		rhs_index.append(end)
		beg = end
	out = {
	'mrca' : mrca,
	'time' : time,
	'lhs.position' : lhs_position,
	'rhs.position' : rhs_position,
	'lhs.index' : lhs_index,
	'rhs.index' : rhs_index }
	return out
