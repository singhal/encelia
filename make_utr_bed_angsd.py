import re

seq_file = '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta'

f = open(seq_file, 'r')

for l in f:
	l = l.rstrip()
	if re.search('>', l):
		seq = f.next().rstrip()
		seqlen = len(seq)
		c = re.search('>(\S+)', l).group(1)
		if re.search('5u(\d+)', l):
			end = re.search('5u(\d+)', l).group(1)
			print '%s:1-%s' % (c, end)
		if re.search('3u',l):
			start = re.search('3u(\d+)', l).group(1)
			print '%s:%s-%s' % (c,start,seqlen)
f.close()

