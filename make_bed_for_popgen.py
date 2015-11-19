import re
from itertools import groupby
import gzip

assembly = '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta'
depth1 = '/home/ssinghal/encelia/variants/depth/Encelia_palmeri.raw_depth.ngm.ancestral_ref.out.gz'
depth2 = '/home/ssinghal/encelia/variants/depth/Encelia_ventorum.raw_depth.ngm.ancestral_ref.out.gz'
# depth1 = 'test1.out'
# depth2 = 'test2.out'
# the number of individuals who need to have that depth
min_depth = 5
# this is out of 5 ind 
min_ind = 4
outfile = '/home/ssinghal/encelia/variants/depth/popgenomics.bed'

###############
# subroutines #
###############

def get_sites(sites, depth, id, min_depth, min_ind):
	f = gzip.open(depth, 'r')
	for l in f:
		d = re.split('\t', l.rstrip())
		if d[0] in sites:
			num_ind = len(filter(lambda x: int(x) > min_depth, d[2:]))
			if num_ind >= min_ind:
				sites[d[0]][id][int(d[1])] = 1
	f.close()
	return sites

def sub(x):
    return x[1] - x[0]

# http://stackoverflow.com/questions/9470611/how-to-do-an-inverse-range-i-e-create-a-compact-range-based-on-a-set-of-numb/9471386#9471386
def get_ranges(a):
	ranges = []
	for k, iterable in groupby(enumerate(sorted(a)), sub):
	     rng = list(iterable)
	     if len(rng) == 1:
	         s = [rng[0][1] - 1, rng[0][1]]
	     else:
	         s = [rng[0][1] - 1, rng[-1][1]]
	     ranges.append(s)
	return ranges

################
# run the code #
################

sites = {}
a = open(assembly, 'r')
for l in a:
	if re.search('gs', l):
		contig = re.search('>(\S+).*gs', l).group(1)
		sites[contig] = {'1': {}, '2': {}}
a.close()

sites = get_sites(sites, depth1, '1', min_depth, min_ind)
sites = get_sites(sites, depth2, '2', min_depth, min_ind)		

out = open(outfile, 'w')
for c in sites:
	shared_sites = set(sites[c]['1'].keys()).intersection(sites[c]['2'].keys())
	shared_sites = list(shared_sites)
	if shared_sites:
		for (i, j) in get_ranges(shared_sites):
			out.write('%s\t%s\t%s\n' % (c, i, j))
out.close()
