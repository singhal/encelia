import re
import pandas as pd
import glob

# we will filter bed file to remove contigs with excess heterozygosity

hetero_file = '/home/ssinghal/encelia/analysis/heterozygosity/inbreeding_fs.palmeri_ventorum.ancestral_ref.csv'
bed_file = '/home/ssinghal/encelia/variants/depth/popgenomics.bed'
out = '/home/ssinghal/encelia/variants/depth/popgenomics.filtered_fs.bed'

d = pd.read_csv(hetero_file)
d['ditch'] = 'keep'
d.ix[(d.avg_fs_ventorum < -0.4), 'ditch'] = 'ditch'
d.ix[(d.avg_fs_palmeri < -0.4), 'ditch'] = 'ditch'
ditch_contigs = d[d.ditch == 'ditch'].contig.tolist()

f = open(bed_file, 'r')
o = open(out, 'w')

count = 0

for l in f:
	d = re.split('\s+', l.rstrip())
	if d[0] not in ditch_contigs:
		o.write(l)
		count += (int(d[2]) - int(d[1]))

f.close()
o.close()
print count
