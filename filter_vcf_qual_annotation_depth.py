import re
import glob
import gzip
import subprocess

vcf_files = glob.glob('/home/ssinghal/encelia/variants/vcfs/*ngm*')
seq = '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta'
bed = '/home/ssinghal/encelia/variants/depth/popgenomics.bed'

anno = {}
s = open(seq, 'r')
for l in s:
	if re.search('gs', l):
		contig = re.search('>(\S+)', l).group(1)
		anno[contig] = 1
s.close()

for vcf_file in vcf_files:
	out_vcf = vcf_file.replace('raw.vcf.gz', 'annotated.filtered_qual.vcf')

	o = open(out_vcf, 'w')
	v = gzip.open(vcf_file, 'r')

	for l in v:
		l = l.rstrip()
		d = re.split('\t', l)
		if re.search('^#', l):
			o.write(l + '\n')
		else:
			if not re.search('INDEL', l):
				ditch = False
				if d[0] not in anno:
					ditch = True
				if re.search('MQB', l):
					qual = float(re.search('MQB=([\d|\-|\.]+)', l).group(1))
					if qual <= 0:
						ditch = True
				if re.search('RPB', l):
					qual = float(re.search('RPB=([\d|\-|\.]+)', l).group(1))
					if qual <= 0:
						ditch = True
				if d[5]:
					qual = float(d[5])
					if qual <= 10:
						ditch = True
				if ditch == False:
					o.write(l + '\n')
	o.close()
	v.close()

	out_vcf2 = out_vcf.replace('qual.vcf', 'qual.high_cov.vcf')
	subprocess.call('~/bin/bedtools2/bin/bedtools intersect -a %s -b %s -sorted > %s' % (out_vcf, bed, out_vcf2), shell=True)
