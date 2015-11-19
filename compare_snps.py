import re

vcf1 = '/home/ssinghal/encelia/bam_files/Encelia_ventorum.bwa.samtools.vcf'
vcf2= '/home/ssinghal/encelia/bam_files/Encelia_ventorum.bowtie.samtools.vcf'

snps = {}
f = open(vcf1, 'r')
for l in f:
	if not re.search('^#', l):
		d = re.split('\t', l)
		if d[0] not in snps:
			snps[d[0]] = {}
		snps[d[0]][d[1]] = 'in1'
f.close()

f = open(vcf2, 'r')
for l in f:
	if not re.search('^#', l):
                d = re.split('\t', l)
                if d[0] not in snps:
                        snps[d[0]] = {}
		if d[1] in snps[d[0]]:
			snps[d[0]][d[1]] = 'both'
		else:
			snps[d[0]][d[1]] = 'in1'
f.close()

o = open('compare_snps.csv', 'w')
o.write('contig,in_both,only_in_bwa,only_in_bowtie\n')
for c in snps:
	num_both = 0
	num_1 = 0
	num_2 = 0
	
	for pos in snps[c]:
		if snps[c][pos] == 'in1':
			num_1 += 1
		elif snps[c][pos] == 'in2':
			num_2 += 1
		else:
			num_both += 1
	o.write('%s,%s,%s,%s\n' % (c, num_both, num_1, num_2))
o.close()

