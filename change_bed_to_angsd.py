import re

d = '/home/ssinghal/encelia/analysis/popgenomics.bed'
out = '/home/ssinghal/encelia/analysis/coverage.angsd_regions.txt'

f = open(d, 'r')
o = open(out, 'w')

for l in f:
	d = re.split('\t', l.rstrip())
	o.write('%s:%s-%s\n' % (d[0], int(d[1]) - 1, d[2]))

f.close()
o.close()
