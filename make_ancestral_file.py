import re
import gzip
import sys
import operator

seqfile = '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta'
aa_file = '/home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.mpileup.out.gz'
#seqfile = '/home/ssinghal/encelia/annotation/test.fasta'
#aa_file = '/home/ssinghal/encelia/variants/ancestral_allele/test.gz'
cutoff = 2

seq = {}
seqf = open(seqfile, 'r')
for l in seqf:
	l = l.rstrip()
	if re.search('>', l):
		c = re.search('>(\S+)', l).group(1)
		seq[c] = list(seqf.next().lower().rstrip())
seqf.close()

def parse_align(bp, ref):
	bp = bp.upper()

	# need to remove indels because that can lead to artificial counts
	nums = [int(x) for x in re.findall('[\-|\+](\d+)', bp)]
	for num in nums:
		bp = re.sub('[\+|\-]%s[A|T|C|G]{%s}' % (num, num), '', bp)
		
	counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
	
	for base in counts.keys():
		counts[base] = bp.count(base) 
	counts[ref] = bp.count('.') + bp.count(',')

	if sum(counts.values()) < cutoff:
		return 'N'
	else:
		countsum = float(sum(counts.values()))
		for (base, count) in counts.items():
			counts[base] = count / countsum
		sorted_counts = sorted(counts.items(), key=operator.itemgetter(1))
		if sorted_counts[3][1] > 0.5:
			return sorted_counts[3][0]
		else:
			return 'N'


f = gzip.open(aa_file)
for ix, l in enumerate(f):
	l = l.rstrip()
	d = re.split('\t', l)
	if d[0] in seq:
		# these are the alignment info
		bases = [parse_align(bp, d[2]) for bp in d[4:][::3]]
		bases = filter(lambda x: x != 'N', bases)

		aa = None
	
		if bases:
			if len(set(bases)) == 1:
				aa = bases[0]
			else:
				counts = {}
				for base in set(bases):
					counts[base] = bases.count(base)
				sorted_counts = sorted(counts.items(), key=operator.itemgetter(1))
				if sorted_counts[-1][1] >= 2:
					aa = sorted_counts[-1][0]
		if aa:
			seq[d[0]][int(d[1]) - 1] = aa
f.close()

o = open('/home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa', 'w')
for c in seq:
	o.write('>%s\n%s\n' % (c, ''.join(seq[c])))
o.close()
