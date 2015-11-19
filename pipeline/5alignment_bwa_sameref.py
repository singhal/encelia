import re
import subprocess
import os
import operator


specieses = {	'Encelia_ventorum': {	'ind': ['CD011', 'CD012', 'CD013', 'CD014', 'CD015'],
					'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' }, 
		'Encelia_palmeri': {	'ind': ['CD659', 'CD662', 'CD663', 'CD664', 'CD661'],
					'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
	}				

'''
specieses = {   'Encelia_asperifolia': {'ind': ['CD104'],
                                        'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
                'Encelia_canescens': {  'ind': ['CD101'],
                                        'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
                'Encelia_farinosa': {   'ind': ['CD103'],
                                        'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
                'Encelia_frutescens': { 'ind': ['CD102'],
                                        'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' }
	}
'''

np = 16
bwa = '/home/ssinghal/bin/bwa'
samtools = '/home/ssinghal/bin/samtools/samtools'
bcftools_dir = '/home/ssinghal/bin/bcftools/'
map_dir = '/home/ssinghal/encelia/variants/bam_files/'
platypus = 'python /home/ssinghal/bin/Platypus_0.7.9.1/Platypus.py'
gatk = '/home/ssinghal/bin/GenomeAnalysisTK.jar'
picard_dir = '/home/ssinghal/bin/picard-tools/'

def make_bwa(seq, bwa):
	if not os.path.isfile(seq + '.bwt'):
		subprocess.call('%s index %s' % (bwa, seq), shell=True)
	return seq 


def run_bwa(map_dir, ind, seq, read_1, read_2, read_u, bwa, np, samtools):
	final = '%s%s.bwa_sameref.bam' % (map_dir, ind)
	
	out1 = '%s%s1bwa' % (map_dir, ind)
	out2 = '%s%s2bwa' % (map_dir, ind)
	tmp = '%s%s.tmp_bwa.bam' % (map_dir, ind)
	finalout = final.replace('.bam', '')
		
	if not os.path.isfile(final):
		subprocess.call("%s mem -t %s %s %s %s > %s.sam" % (bwa, np, seq, read_1, read_2, out1), shell=True)
		subprocess.call("%s mem -t %s %s %s > %s.sam" % (bwa, np, seq, read_u, out2), shell=True)
		subprocess.call("%s view -bS %s.sam > %s.bam" % (samtools, out1, out1), shell=True)
		subprocess.call("%s view -bS %s.sam > %s.bam" % (samtools, out2, out2), shell=True)
		subprocess.call("%s merge %s %s.bam %s.bam" % (samtools, tmp, out1, out2), shell=True)
		subprocess.call("%s sort %s %s" % (samtools, tmp, finalout), shell=True)
		subprocess.call("rm %s.sam %s.sam %s.gz %s.bam %s.bam %s" % (out1, out2, out1, out1, out2, tmp), shell=True)
	return final

def run_samtools(seq, bams, samtools, bcftools_dir, map_dir, species):
	out_bcf = '%s%s.bwa_sameref.bcf' % (map_dir, species)
	out_vcf = '%s%s.bwa_sameref.samtools.vcf' % (map_dir, species)
	if not os.path.isfile(out_vcf):
		subprocess.call('%s mpileup -B -uf %s %s | %sbcftools call -m -O b -v - > %s' % (samtools, seq, ' '.join(bams), bcftools_dir, out_bcf), shell=True)
		subprocess.call('%sbcftools view %s | %svcfutils.pl varFilter > %s' % (bcftools_dir, out_bcf, bcftools_dir, out_vcf), shell=True)
	return out_vcf

def get_fixed(species_vcf):
	var = {}
	f = open(species_vcf, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			alleles = [d[3]] + re.split(',', d[4])
			indel = False
			for allele in alleles:
				if len(allele) > 1:
					indel = True		
			
			allele_counts = {}
			for ix, allele in enumerate(alleles):
				allele_counts[ix] = len(re.findall('%s/' % ix, l)) + len(re.findall('/%s' % ix, l))
			missing = len(re.findall('\./', l)) + len(re.findall('/\.', l))		
			
			# don't want missing data to be too high
			if missing / float(sum(allele_counts.values())) < 0.5:
				# if alternate allele is more than 80%
				for ix in allele_counts:
					allele_counts[ix] = allele_counts[ix] / float(sum(allele_counts.values()))
				high_allele = max(allele_counts.iteritems(), key=operator.itemgetter(1))[0]
				if high_allele != 0 and allele_counts[high_allele] >= 0.8:
					if d[0] not in var:
						var[d[0]] = dict()
					var[d[0]][int(d[1])] = {'ref': d[3], 'alt': alleles[int(high_allele)], 'indel': indel}
	f.close()
	return var
					
def get_seq(seq):
	f = open(seq, 'r')
	contig = ''
	seqs = {}
	for l in f:
		if re.search('>', l):
			contig = re.search('>(\S+)', l).group(1)
			seqs[contig] = ''
		else:
			seqs[contig] += l.rstrip()
	f.close()
	for contig, seq in seqs.items():
		seqs[contig] = list(seq)
	return seqs

			
def mutate_seqs(out, seq, fixed_var):
	o = open(out, 'w')
	for contig in seq:
		if contig in fixed_var:
			for pos in sorted(fixed_var[contig]):
				seq[contig][pos-1] = fixed_var[contig][pos]['alt']
		o.write('>%s\n%s\n' % (contig, ''.join(seq[contig])))
	o.close()
	return
		
def run_platypus(platypus, seq, bams, samtools, map_dir, species, np):
        out = '%s%s.bwa.platypus.vcf' % (map_dir, species)
        if not os.path.isfile(out):
		if not os.path.isfile(seq + '.fai'):                
			subprocess.call("%s faidx %s" % (samtools, seq), shell=True)
                for bam in bams:
			if not os.path.isfile(bam + '.bai'):               
				subprocess.call("%s index %s" % (samtools, bam), shell=True)
                	subprocess.call("%s callVariants --bamFiles=%s --refFile=%s --output=%s --nCPU %s --maxVariants 20" % (platypus, ','.join(bams), seq, out, np), shell=True)


for species in specieses:
	bams = []
	for ind in specieses[species]['ind']:
		seq = specieses[species]['seq']
		read_1 = '/home/ssinghal/encelia/trimmed_reads/%s/%s/%s_1_final.fastq' % (species, ind, ind)
		read_2 = read_1.replace('_1_', '_2_')
		read_u = read_1.replace('_1_', '_u_')

		seqout = make_bwa(seq, bwa)
		bam = run_bwa(map_dir, ind, seqout, read_1, read_2, read_u, bwa, np, samtools)
		bams.append(bam)
	species_vcf = run_samtools(seq, bams, samtools, bcftools_dir, map_dir, species)
	# run_platypus(platypus, seq, bams, samtools, map_dir, species, np)
	# fixed_var = get_fixed(species_vcf)
	# sequences = get_seq(seq)
	# out = seq + '.fixed'
	# mutate_seqs(out, sequences, fixed_var)
