import re

vcf_file1 = '/home/ssinghal/encelia/variants/Encelia_palmeri.bwa_sameref.samtools.annotated.filteredSNPs.vcf'
vcf_file2 = '/home/ssinghal/encelia/variants/Encelia_ventorum.bwa_sameref.samtools.annotated.filteredSNPs.vcf'

def parse_vcf(file):
	f = open(file, 'r')
	snps = {}
	for l in f:
		if not re.search('^#', l):
			genos = re.findall('(\S)\/', l) + re.findall('\/(\S)', l)
			genos = filter(lambda x: x != '.', genos)
			ref = float(genos.count('0'))/ float(len(genos))
			
			d = re.split('\t', l)
			if d[0] not in snps:
				snps[d[0]] = {}
			alt = re.split(',', d[4])
			if not isinstance(alt, list):
				alt = [alt]
			snps[d[0]][d[1]] = {'alt': alt, 'freq': 1 - ref}
	f.close()
	return snps	
				

snps1 = parse_vcf(vcf_file1)
snps2 = parse_vcf(vcf_file2)

print 'contig,pos,type,freq'

for c in snps1:
	if c in snps2:
		for pos in snps1[c]:
			if pos in snps2[c]:
				if snps1[c][pos]['alt'] == snps2[c][pos]['alt']:
					freq = abs(snps1[c][pos]['freq'] - snps2[c][pos]['freq'])
					print '%s,%s,shared,%s' % (c, pos, freq)
				else:
					print '%s,%s,palmeri,%s' % (c, pos, snps1[c][pos]['freq'])
					print '%s,%s,ventorum,%s' % (c, pos, snps2[c][pos]['freq'])
			else:
				print '%s,%s,palmeri,%s' % (c, pos, snps1[c][pos]['freq'])
	else:
		for pos in snps1[c]:
			print '%s,%s,palmeri,%s' % (c, pos, snps1[c][pos]['freq'])


for c in snps2:
        if c in snps1:
                for pos in snps2[c]:
                        if pos not in snps1[c]:
				print '%s,%s,ventorum,%s' % (c, pos, snps2[c][pos]['freq'])
        else:
		for pos in snps2[c]:
	                print '%s,%s,ventorum,%s' % (c, pos, snps2[c][pos]['freq'])
