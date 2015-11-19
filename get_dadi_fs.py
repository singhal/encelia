import re
import random
import numpy as np
import copy
import subprocess 

seqfile = '/home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa'
# vcf1 = '/home/ssinghal/encelia/variants/test_palmeri.vcf'
# vcf2 = '/home/ssinghal/encelia/variants/test_ventorum.vcf'

vcf1 = '/home/ssinghal/encelia/variants/vcfs/palmeri.ngm.annotated.filtered_qual.high_cov.filtered_fs.vcf'
vcf2 = '/home/ssinghal/encelia/variants/vcfs/ventorum.ngm.annotated.filtered_qual.high_cov.filtered_fs.vcf'
num_ind = 5
out = '/home/ssinghal/encelia/analysis/dadi/palmeri_ventorum.same_ref.fs'

def get_variants(vcf):
	var = {}
	f = open(vcf, 'r')
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
				for geno in d[9:]:
					geno = re.search('^(\S\/\S)', geno).group(1)
					geno = re.split('\/', geno)
					genos += geno
				genos = [x for x in genos if not re.search('\.', x)]
					
				counts = {}
				for ix, allele in enumerate(alleles):
					count = genos.count(str(ix))
					if count > 0:
						counts[allele] = count
				if len(counts) > 0:
					if d[0] not in var:
						var[d[0]] = {}
					var[d[0]][int(d[1])] = counts

	return var


def get_sequence(seqfile, contig):
	seqcall = subprocess.Popen('~/bin/samtools/samtools faidx %s %s' % (seqfile, contig), shell=True, stdout=subprocess.PIPE)
	seq = ''
	for l in seqcall.stdout:
		if not re.match('>', l):
			seq += l.rstrip()
	seq = list(seq)
	return seq


def get_values(values, call, allele1, allele2):
	if allele1 not in call and allele2 not in call:
		values['allele1_counts'].append(10)
		values['allele2_counts'].append(0)
	else:
		if allele1 in call:
			values['allele1_counts'].append(call[allele1])
			if allele2 in call:
				values['allele2_counts'].append(call[allele2])
			else:
				values['allele2_counts'].append(10 - call[allele1])
      		else:
			if allele2 in call:
				values['allele1_counts'].append(10 - call[allele2])
				values['allele2_counts'].append(call[allele2])
	return values


def get_allele_counts(o, call1, call2, seq, pos, contig):
	alleles = list(set(call1.keys() + call2.keys()))
	aa_triplet = seq[pos-2:pos+1]

	values = {'allele1': None, 'allele1_counts': [], 'allele2': None, 'allele2_counts': [], 'aa_trip': ''.join(aa_triplet)}

	# don't want to deal with undefined aa or multiallelics
	if aa_triplet[1] not in ['a', 't', 'c', 'g'] or len(alleles) > 2:
		allele1 = aa_triplet[1]
		# don't want to deal with homoplasy
		if len(alleles) == 2 and allele1 in alleles:
			values['allele1'] = allele1
			allele2 = [x for x in alleles if x != allele1][0]
			values['allele2'] = allele2
			values = get_values(values, call1, allele1, allele2)
			values = get_values(values, call2, allele1, allele2)
		# fixed for the other
		elif len(alleles) == 1 and allele1 not in alleles:
			values['allele1'] = allele1
			allele2 = alleles[0]
			values['allele2'] = allele2
			values = get_values(values, call1, allele1, allele2)
                        values = get_values(values, call2, allele1, allele2)
	
	if values['allele1']:
		o.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (values['aa_trip'], values['aa_trip'], values['allele1'], values['allele1_counts'][0], \
				values['allele1_counts'][1], values['allele2'], values['allele2_counts'][0], values['allele2_counts'][1], contig, pos))
	

def make_fs(out, var1, var2, num_ind, seqfile):
	o = open(out, 'w')
	o.write('outgroup1\toutgroup2\tAllele1\tpalmeri\tventorum\tAllele1\tpalmeri\tventorum\tcontig\tposition\n')

	for contig in var1:
		seq = get_sequence(seqfile, contig)

		if contig in var2:
			# could possibly have shared snps
			for pos in var1[contig]:
				if pos in var2[contig]:
					# shared snp
					call1 = var1[contig][pos]
					call2 = var2[contig][pos]
					get_allele_counts(o, call1, call2, seq, pos, contig)
				else:
					call1 = var1[contig][pos]
	                                get_allele_counts(o, call1, {}, seq, pos, contig)

			for pos in var2[contig]:
				if pos not in var1[contig]:
					call2 = var2[contig][pos]
                                	get_allele_counts(o, {}, call2, seq, pos, contig)

		else:
			# all snps will only be in sp1
			for pos in var1[contig]:
				call1 = var1[contig][pos]
				get_allele_counts(o, call1, {}, seq, pos, contig)
				
	for contig in var2:
		if contig not in var1:
			seq = get_sequence(seqfile, contig)
			# these are the unique positions that i haven't considered yet
			for pos in var2[contig]:
				call2 = var2[contig][pos]
                                get_allele_counts(o, {}, call2, seq, pos, contig)
	o.close() 
				
palmeri_var = get_variants(vcf1)
ventorum_var = get_variants(vcf2)
make_fs(out, palmeri_var, ventorum_var, num_ind, seqfile)



