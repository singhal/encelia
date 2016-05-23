import re
import gzip

vcf1 = '/home/ssinghal/encelia/variants/vcfs/palmeri.ngm.annotated.filtered_qual.high_cov.vcf.gz'
vcf2 = '/home/ssinghal/encelia/variants/vcfs/ventorum.ngm.annotated.filtered_qual.high_cov.vcf.gz'
bed = '/home/ssinghal/encelia/variants/depth/popgenomics.bed'
# vcf1 = '/home/ssinghal/encelia/variants/vcfs/test1.vcf.gz'
# vcf2 = '/home/ssinghal/encelia/variants/vcfs/test2.vcf.gz'
nind1 = 10
nind2 = 10

snps1 = {}
snps2 = {}

def get_SNPs(vcf, snps):
	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('#', l):
			d = re.split('\t', l.rstrip())
			alleles = re.split(',', d[4])
			indel = False
			for allele in alleles:
				if len(allele) > 1:
					indel = True
			if not indel:
				orig_genos = [re.search('^(\S/\S)', x).group(1) for x in d[9:]]
				genos = []
				for x in orig_genos:
					genos += re.split('/', x)
				
				if d[0] not in snps:
					snps[d[0]] = {}
				snps[d[0]][d[1]] = genos
	f.close()
	return snps


def contig_snps(snps, c):
	if c in snps:
		return snps[c].keys()
	else:
		return []


def get_genos(snps, c, snp, nind):
	genos = []
	if c in snps:
		if snp in snps[c]:
			genos = snps[c][snp]
	if not genos:
		genos = ['0'] * nind
	return genos


def pi_within(genos):
	pi = 0
	genos = [x for x in genos if x != '.']
	inds = float(len(genos))
	alleles = set(genos)
	for allele1 in alleles:
		for allele2 in alleles:
			diff = 0
			if allele1 != allele2:
				diff = 1
			pi += genos.count(allele1) / inds * genos.count(allele2) / inds * diff
	return pi


def pi_btn(genos1, genos2):
        pi = 0

        genos1 = [x for x in genos1 if x != '.']
	genos2 = [x for x in genos2 if x != '.']

        inds1 = float(len(genos1))
	inds2 = float(len(genos2))

        alleles = set(genos1 + genos2)

        for allele1 in alleles:
                for allele2 in alleles:
                        diff = 0
                        if allele1 != allele2:
                                diff = 1
                        pi += genos1.count(allele1) / inds1 * genos2.count(allele2) / inds2 * diff
        return pi


def get_divergence(snps1, snps2):
	div = {}
	contigs = set(snps1.keys() + snps2.keys())

	for c in contigs:
		div[c] = {'pi': 0, 'pi1': 0, 'pi2': 0}
		snps = set(contig_snps(snps1, c) + contig_snps(snps2, c))
		for snp in snps:
			genos1 = get_genos(snps1, c, snp, nind1)
			genos2 = get_genos(snps2, c, snp, nind2)

			div[c]['pi1'] += pi_within(genos1)
			div[c]['pi2'] += pi_within(genos2)
			div[c]['pi'] += pi_btn(genos1, genos2)

	return div


def get_denom(bed):
	denom = {}
	f = open(bed, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		if d[0] not in denom:
			denom[d[0]] = 0
		denom[d[0]] += int(d[2]) - int(d[1])
	f.close()
	return denom


def report_div(denom, div):
	print 'contig,denom,pi1,pi2,pi'
	keys = ['pi1', 'pi2', 'pi']
	total = {'pi1': 0, 'pi2': 0, 'pi': 0}
	for c in denom:
		if c in div:
			vals = [div[c][key] / float(denom[c]) for key in keys]
			for key in keys:
				total[key] += div[c][key]
		else:
			vals = [0, 0, 0]
		print '%s,%s,%s' % (c, denom[c], ','.join([str(x) for x in vals]))
	total_denom = sum(denom.values())
	vals = [total[key] / float(total_denom) for key in keys]
	
	print 'all,%s,%s' % (total_denom, ','.join([str(x) for x in vals]))

snps1 = get_SNPs(vcf1, snps1)
snps2 = get_SNPs(vcf2, snps2)
div = get_divergence(snps1, snps2)
denom = get_denom(bed)
report_div(denom, div)
