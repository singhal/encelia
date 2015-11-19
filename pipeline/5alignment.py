import re
import subprocess
import os
import operator

specieses = {	'Encelia_ventorum': {	'ind': ['CD011', 'CD012', 'CD013', 'CD014', 'CD015'],
					'seq': '/home/ssinghal/encelia/annotation/ventorum.annotated.fasta' }, 
		'Encelia_palmeri': {	'ind': ['CD659', 'CD662', 'CD663', 'CD664', 'CD661'],
					'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
		'Encelia_asperifolia': {'ind': ['CD104'],
					'seq': '/home/ssinghal/encelia/annotation/asperifolia.annotated.fasta' },
		'Encelia_canescens': {  'ind': ['CD101'],
					'seq': '/home/ssinghal/encelia/annotation/canescens.annotated.fasta' },
		'Encelia_farinosa': {	'ind': ['CD103'],
					'seq': '/home/ssinghal/encelia/annotation/farinosa.annotated.fasta' },
		'Encelia_frutescens': {	'ind': ['CD102'],
					'seq': '/home/ssinghal/encelia/annotation/frutescens.annotated.fasta' }
	}	
			
np = 16
java = '/home/ssinghal/bin/jre1.7.0_67/bin/java'
bowtie = '/home/ssinghal/bin/bowtie2'
bowtie_build = '/home/ssinghal/bin/bowtie2-build'
samtools = '/home/ssinghal/bin/samtools/samtools'
bcftools_dir = '/home/ssinghal/bin/bcftools/'
map_dir = '/home/ssinghal/encelia/bam_files/'
platypus = 'python /home/ssinghal/bin/Platypus_0.7.9.1/Platypus.py'
gatk = '/home/ssinghal/bin/GenomeAnalysisTK.jar'
picard_dir = '/home/ssinghal/bin/picard-tools/'

def make_bowtie(seq, bowtie_build):
	# if not os.path.isfile(seq + '.1.bt2'):
	#	subprocess.call('%s %s %s' % (bowtie_build, seq, seq), shell=True)
	return seq 


def run_bowtie(map_dir, ind, seq, read_1, read_2, read_u, bowtie, np, samtools):
	final = '%s%s.bowtie.bam' % (map_dir, ind)
	
	out1 = '%s%s1' % (map_dir, ind)
	out2 = '%s%s2' % (map_dir, ind)
	tmp = '%s%s.tmp.bam' % (map_dir, ind)
	finalout = final.replace('.bam', '')
		
	if not os.path.isfile(final):
		subprocess.call("%s -p %s -x %s -1 %s -2 %s -S %s.sam -X 600" % (bowtie, np, seq, read_1, read_2, out1), shell=True)
		subprocess.call("%s -p %s -x %s -U %s -S %s.sam" % (bowtie, np, seq, read_u, out2), shell=True)
		subprocess.call("%s view -bS %s.sam > %s.bam" % (samtools, out1, out1), shell=True)
		subprocess.call("%s view -bS %s.sam > %s.bam" % (samtools, out2, out2), shell=True)
		subprocess.call("%s merge %s %s.bam %s.bam" % (samtools, tmp, out1, out2), shell=True)
		subprocess.call("%s sort %s %s" % (samtools, tmp, finalout), shell=True)
		subprocess.call("rm %s.sam %s.sam %s.gz %s.bam %s.bam %s" % (out1, out2, out1, out1, out2, tmp), shell=True)
	return final

def run_samtools(seq, bams, samtools, bcftools_dir, map_dir, species):
	out_bcf = '%s%s.raw.bcf' % (map_dir, species)
	out_vcf = '%s%s.bowtie.samtools.vcf' % (map_dir, species)
	if not os.path.isfile(out_vcf):
		subprocess.call('%s mpileup -C 50 -uf %s %s | %sbcftools call -m -O b -v - > %s' % (samtools, seq, ' '.join(bams), bcftools_dir, out_bcf), shell=True)
		subprocess.call('%sbcftools view %s | %svcfutils.pl varFilter -w 0 > %s' % (bcftools_dir, out_bcf, bcftools_dir, out_vcf), shell=True)
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
	out = '%s%s.bowtie.platypus.vcf' % (map_dir, species)
	if not os.path.isfile(out):
		subprocess.call("%s faidx %s" % (samtools, seq), shell=True)
		for bam in bams:
			subprocess.call("%s index %s" % (samtools, bam), shell=True)
		subprocess.call("%s callVariants --bamFiles=%s --refFile=%s --output=%s --nCPU %s --maxVariants 20" % (platypus, ','.join(bams), seq, out, np), shell=True)


def run_gatk(java, gatk, picard_dir, seq, bams, samtools, map_dir, species, np):
	mid_bams = []
	for bam in bams:
		lib = re.search('(CD\d+)', bam).group(1)
		out_bam1 = '%s%s.rg.bam' % (map_dir, lib)
		out_bam2 = '%s%s.rg.dup.bam' % (map_dir, lib)
		out_bam3 = '%s%s.rg.dup.realigned.bam' % (map_dir, lib)
		out_bam4 = '%s%s.rg.dup.realigned.mateFixed.bam' % (map_dir, lib)

		if not os.path.isfile(out_bam1):
			subprocess.call("%s -Xmx6g -jar %sAddOrReplaceReadGroups.jar I=%s O=%s SO=coordinate RGID=Encelia RGLB=%s RGPL=Illumina RGPU=GAIIx RGSM=%s MAX_RECORDS_IN_RAM=1000000 TMP_DIR=/home/ssinghal/tmp/" % (java, picard_dir, bam, out_bam1, lib, lib), shell=True)
		if not os.path.isfile(out_bam2):
			subprocess.call("%s -Xmx6g -jar %sMarkDuplicates.jar I=%s O=%s  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=%s.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000" % (java, picard_dir, out_bam1, out_bam2, lib), shell=True)
		if not os.path.isfile(out_bam3):
			subprocess.call("%s -Xmx2g -jar %s -T RealignerTargetCreator -R %s -I %s -o %s.forIndelRealigner.intervals" % (java, gatk, seq, out_bam2, lib), shell=True)
			subprocess.call("%s -Xmx4g -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s.forIndelRealigner.intervals -o %s" % (java, gatk, seq, out_bam2, lib, out_bam3), shell=True)
		if not os.path.isfile(out_bam4):
			subprocess.call("%s -Xmx4g -jar %sFixMateInformation.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/ssinghal/tmp/" % (java, picard_dir, out_bam3, out_bam4), shell=True)
		mid_bams.append(out_bam4)

	raw_vcf = '%s%s.bowtie.gatk.raw.vcf' % (map_dir, species)
	mid_bam_str = '-I ' + ' -I '.join(mid_bams)
	subprocess.call("%s -Xmx4g -jar %s -T UnifiedGenotyper -R %s %s -glm BOTH -minIndelCnt 1 -mbq 20 -hets 0.035 -out_mode EMIT_ALL_SITES -o %s -nct %s" % (java, gatk, seq, mid_bam_str, raw_vcf, np), shell=True)
	
	final_bams = []
	for bam in mid_bams:
		out_bam5 = '%s%s.rg.dup.realigned.mateFixed.recal.bam' % (map_dir, lib)
		subprocess.call("%s -Xmx4g -jar %s -T BaseRecalibrator -R %s %s -knownSites %s -o %s.recal.grp -nct %s" % (java, gatk, seq, mid_bam_str, raw_vcf, lib, np), shell=True)
		subprocess.call("%s -Xmx4g -jar %s -T PrintReads -R %s %s -BQSR %s.recal.grp -o %s -nct %s" % (java, gatk, seq, mid_bam_str, lib, out_bam5, np), shell=True)
		final_bams.append(out_bam5)
	
	final_bam_str = '-I ' + ' -I '.join(final_bams)
	hc_vcf = '%s%s.bowtie.gatk.raw.bqsr.vcf' % (map_dir, species)
	final_vcf = '%s%s.bowtie.gatk.vcf' % (map_dir, species)	

	subprocess.call("%s -Xmx8g -jar %s -T HaplotypeCaller -R %s %s -stand_call_conf 20.0 -stand_emit_conf 20.0 -o %s -hets 0.035 -mbq 20 --recoverDanglingHeads" % (java, gatk, seq, final_bam_str, hc_vcf), shell=True)
	subprocess.call('%s -Xmx4g -jar %s -T VariantFiltration -R %s -V %s -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o %s' % (java, gatk, seq, hc_vcf, final_vcf), shell=True)


for species in specieses:
	bams = []
	for ind in specieses[species]['ind']:
		seq = specieses[species]['seq']
		read_1 = '/home/ssinghal/encelia/trimmed_reads/%s/%s/%s_1_final.fastq' % (species, ind, ind)
		read_2 = read_1.replace('_1_', '_2_')
		read_u = read_1.replace('_1_', '_u_')

		seqout = make_bowtie(seq, bowtie_build)
		bam = run_bowtie(map_dir, ind, seqout, read_1, read_2, read_u, bowtie, np, samtools)
		bams.append(bam)
	# species_vcf = run_samtools(seq, bams, samtools, bcftools_dir, map_dir, species)
	# run_platypus(platypus, seq, bams, samtools, map_dir, species, np)
	run_gatk(java, gatk, picard_dir, seq, bams, samtools, map_dir, species, np)

	# fixed_var = get_fixed(species_vcf)
	# sequences = get_seq(seq)
	# out = seq + '.fixed'
	# mutate_seqs(out, sequences, fixed_var)
