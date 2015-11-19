import re
import pandas as pd

# we will filter dadi file both by contigs with excess heterozygosity
# and snps with missing data

original_file = '/home/ssinghal/encelia/analysis/popgen/popgenomics.bed'
hetero_file = '/home/ssinghal/encelia/analysis/dadi/palmeri_ventorum.same_ref.fs_inbreeding.csv'

d = pd.read_csv(hetero_file)
d['ditch'] = 'keep'
d.ix[(d.avg_fs_ventorum < -0.4), 'ditch'] = 'ditch'
d.ix[(d.avg_fs_palmeri < -0.4), 'ditch'] = 'ditch'
ditch_contigs = d[d.ditch == 'ditch'].contig.tolist()

f = open(original_file, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	if d[0] not in ditch_contigs:
		print l.rstrip()
