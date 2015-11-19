import re
import numpy as np

vcf1 = '/home/ssinghal/encelia/variants/vcfs/palmeri.ngm.annotated.filtered_qual.high_cov.vcf'
vcf2 = '/home/ssinghal/encelia/variants/vcfs/ventorum.ngm.annotated.filtered_qual.high_cov.vcf'

def get_hetero(vcf):
	f = open(vcf, 'r')
	hetero = {}

	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			
			alleles = [d[3]] + re.split(',', d[4])
                        indel = False
                        for allele in alleles:
                                if len(allele) > 1:
                                        indel = True

                        if not indel:
                                genos = []
				obs_het = 0
                                for geno in d[9:]:
                                        geno = re.search('^(\S\/\S)', geno).group(1)
                                        geno = re.split('\/', geno)
					if geno[0] != geno[1]:
						obs_het += 1
                                        genos += geno
                                genos = [x for x in genos if not re.search('\.', x)]
				obs_het = obs_het / float(len(genos)/2.0)                                        

                                counts = {}
                                for ix, allele in enumerate(alleles):
                                        count = genos.count(str(ix))
                                        if count > 0:
                                                counts[allele] = count

				# only calc heterozygosity for biallelic sites, could calc for more, but not worth it
				if len(counts) == 2:
					count_alleles = counts.values()
					p = count_alleles[0] / float(len(genos))
					q = count_alleles[1] / float(len(genos))
					exp_het = 2 * p * q
					
					fs = (exp_het - obs_het) / exp_het

					if d[0] not in hetero:
						hetero[d[0]] = {'fs': 0, 'num': 0}
					hetero[d[0]]['fs'] += fs
					hetero[d[0]]['num'] += 1

	return hetero


def print_results(hetero1, hetero2):
	print 'contig,avg_fs_palmeri,num_palmeri,avg_fs_ventorum,num_ventorum'
	contigs = list(set(hetero1.keys() + hetero2.keys()))
	
	for c in contigs:
		avg_fs_palmeri = 'nan'
		num_palmeri = 'nan'
		avg_fs_ventorum = 'nan'
		num_ventorum = 'nan'

		if c in hetero1:
			avg_fs_palmeri = hetero1[c]['fs'] / float(hetero1[c]['num'])
			num_palmeri = hetero1[c]['num']
		if c in hetero2:
			avg_fs_ventorum = hetero2[c]['fs'] / float(hetero2[c]['num'])
                        num_ventorum = hetero2[c]['num']	

		print '%s,%s,%s,%s,%s' % (c, avg_fs_palmeri, num_palmeri, avg_fs_ventorum, num_ventorum)

hetero1 = get_hetero(vcf1)
hetero2 = get_hetero(vcf2)
print_results(hetero1, hetero2)
