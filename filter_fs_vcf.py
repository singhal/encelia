import re
import pandas as pd
import glob

# we will filter dadi file both by contigs with excess heterozygosity
# and snps with missing data

hetero_file = '/home/ssinghal/encelia/analysis/heterozygosity/inbreeding_fs.palmeri_ventorum.ancestral_ref.csv'
vcf_files = glob.glob('/home/ssinghal/encelia/variants/vcfs/*high_cov*')

d = pd.read_csv(hetero_file)
d['ditch'] = 'keep'
d.ix[(d.avg_fs_ventorum < -0.4), 'ditch'] = 'ditch'
d.ix[(d.avg_fs_palmeri < -0.4), 'ditch'] = 'ditch'
ditch_contigs = d[d.ditch == 'ditch'].contig.tolist()

for file in vcf_files:
	out = file.replace('.vcf', '.filtered_fs.vcf')
	o = open(out, 'w')
	f = open(file, 'r')
	for l in f:
		d = re.split('\t', l.rstrip())
		if d[0] not in ditch_contigs:
			o.write(l.rstrip() + '\n')
	o.close()
	f.close()
