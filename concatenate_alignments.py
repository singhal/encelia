import glob
import re
import random

files = glob.glob('/home/ssinghal/encelia/analysis/gene_trees/alignments/*aln')

num_aligns = 1

for ix in range(0,num_aligns):
	out = '/home/ssinghal/encelia/analysis/gene_trees/sampled_genes.%s.concatenated_woutgroup.phy' % ix
	random.shuffle(files)
	species = {}
	species['helianthus'] = ''

	for i, file in enumerate(files):
		print i
		if len(species.values()[0]) < 5000000:
			id = ''
			f = open(file, 'r')
			for l in f:
				if re.search('>', l):
					id = re.search('>(\S+)', l).group(1)
					if id not in species:
						species[id] = ''
				else:
					species[id] += l.rstrip()
		else:
			break

	o = open(out, 'w')
	o.write(' %s %s\n' % (len(species), len(species.values()[0])))
	for id in species:
		id2 = id + (20 - len(id)) * ' '
		o.write('%s%s\n' % (id2, species[id]))
	o.close()
	
