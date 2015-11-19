import re
import subprocess
import os
import operator

specieses = {   'Encelia_ventorum': {   'ind': ['CD011', 'CD012', 'CD013', 'CD014', 'CD015'],
                                        'seq': '/home/ssinghal/encelia/annotation/ventorum.annotated.fasta' }, 
                'Encelia_palmeri': {    'ind': ['CD659', 'CD662', 'CD663', 'CD664', 'CD661'],
                                        'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' }
#                'Encelia_asperifolia': {'ind': ['CD104'],
#                                        'seq': '/home/ssinghal/encelia/annotation/asperifolia.annotated.fasta' },
#                'Encelia_canescens': {  'ind': ['CD101'],
#                                        'seq': '/home/ssinghal/encelia/annotation/canescens.annotated.fasta' },
#                'Encelia_farinosa': {   'ind': ['CD103'],
#                                        'seq': '/home/ssinghal/encelia/annotation/farinosa.annotated.fasta' },
#               'Encelia_frutescens': { 'ind': ['CD102'],
#                                        'seq': '/home/ssinghal/encelia/annotation/frutescens.annotated.fasta' }
        }

dir = '/home/ssinghal/encelia/variants/bam_files/'
samtools = '/home/ssinghal/bin/samtools/samtools'


for species in specieses:
	bam_files = ['%s%s.ngm.palmeri_ref.bam' % (dir, ind) for ind in specieses[species]['ind']]
	tmp_out = '%s%s.raw_depth.ngm.palmeri_ref.out' % (dir, species)
	out = '%s%s.depth.ngm.palmeri_ref.csv' % (dir, species) 

	subprocess.call("%s depth %s > %s" % (samtools, ' '.join(bam_files), tmp_out), shell=True)

	'''
	contig = {}
	f = open(tmp_out, 'r')
	for l in f:
		l = l.rstrip()
		d = re.split('\t', l)
		
		if d[0] in contig:
			for depth in d[2:]:
				contig[d[0]]['depth'] += int(depth)
			contig[d[0]]['bp'] += 1
		else:
			contig[d[0]] = {'depth': 0, 'bp': 0}
	f.close()

	o = open(out, 'w')
	o.write('contig,depth\n')
	for c in contig:
		avg_depth = contig[c]['depth'] / float( contig[c]['bp'] * len(bam_files) )
		o.write('%s,%.1f\n' % (c, avg_depth))
	o.close()

	# os.remove(tmp_out)
	'''
