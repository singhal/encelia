import re
import random
import numpy as np
import copy

seqfile = '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta'
# seqfile = '/home/ssinghal/encelia/annotation/test.annotated.fasta'
vcf1 = '/home/ssinghal/encelia/variants/Encelia_palmeri.bwa_sameref.samtools.annotated.filteredSNPs.vcf'
vcf2 = '/home/ssinghal/encelia/variants/Encelia_ventorum.bwa_sameref.samtools.annotated.filteredSNPs.vcf'
# vcf1 = '/home/ssinghal/encelia/variants/test_palmeri.vcf'
# vcf2 = '/home/ssinghal/encelia/variants/test_ventorum.vcf'
num_ind = 5

def get_seq(seqfile):
	seq = {}
	f = open(seqfile, 'r')
	for l in f:
		if re.search('gs', l):
			c = re.search('>(\S+)', l).group(1)
			seq[c] = list(f.next().rstrip())
	f.close()
	return seq

def get_haplo(seq, vcf, num_ind):
	haplo = {}
	for c in seq:
		haplo[c] = {}
		for i in range(num_ind * 2):
			haplo[c][i] = copy.copy(seq[c])
	
	v = open(vcf, 'r')
	for l in v:
		l = l.rstrip()
		if not re.search('^#', l):
			d = re.split('\t', l)
			pos = int(d[1]) - 1
			
			alleles = {}
			for i, allele in enumerate([d[3]] + re.split(',', d[4])):
				alleles[str(i)] = allele
			alleles['.'] = d[3]

			for ix, (geno1, geno2) in enumerate(re.findall('(\S)\/(\S)', l)):
				genos = [geno1, geno2]
				random.shuffle(genos)
		
				ind1 = 2 * ix
				ind2 = 2 * ix + 1
	
				haplo[d[0]][ind1][pos] = alleles[genos[0]]
				haplo[d[0]][ind2][pos] = alleles[genos[1]]
	v.close()
	return haplo

def get_dx(c, num_ind, haplo):
	diffs = []
	for i in range(num_ind * 2):
		for j in range(i+1, num_ind *2):
			num_diff = 0
			for bp1, bp2 in zip(haplo[c][i], haplo[c][j]):
				if bp1 != bp2:
					num_diff += 1
			diffs.append(num_diff)
	seq_diff =  np.mean(diffs) / float(len(haplo[c][0]))
	return seq_diff

def get_dxy(c, num_ind, haplo1, haplo2):
        diffs = []
        for i in range(num_ind * 2):
                for j in range(num_ind * 2):
                        num_diff = 0
                        for bp1, bp2 in zip(haplo1[c][i], haplo2[c][j]):
                                if bp1 != bp2:
                                        num_diff += 1
                        diffs.append(num_diff)
        seq_diff =  np.mean(diffs) / float(len(haplo1[c][0]))
	return seq_diff


seq = get_seq(seqfile)
haplo1 = get_haplo(seq, vcf1, num_ind)
haplo2 = get_haplo(seq, vcf2, num_ind)

out = open('divergence_palmeri_ventorum.csv','w')
out.write('c,dx_palmeri,dx_ventorum,dxy,da\n')
for c in seq:
	dx1 = get_dx(c, num_ind, haplo1)
	dx2 = get_dx(c, num_ind, haplo2)
	dxy = get_dxy(c, num_ind, haplo1, haplo2)	
	
	da = dxy - (dx1 + dx2) / 2.0		
			
	out.write('%s,%.3f,%.3f,%.3f,%.3f\n' % (c, dx1, dx2, dxy, da))
out.close()
	

