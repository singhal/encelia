import re
import subprocess
import os
import operator

specieses = {	'Encelia_ventorum': {	'div': 0.03,
					'ind': ['CD011', 'CD012', 'CD013', 'CD014', 'CD015'],
					'seq': '/home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa' }, 
		'Encelia_palmeri': {	'div': 0.005,
					'ind': ['CD659', 'CD662', 'CD663', 'CD664', 'CD661'],
					'seq': '/home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa' }
		#'Encelia_asperifolia': { 'div': 0.03,
		#			'ind': ['CD104'],
		#			'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
		#'Encelia_canescens': {  'div': 0.03,
		#			'ind': ['CD101'],
		#			'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
		#'Encelia_farinosa': {	'div': 0.03,
		#			'ind': ['CD103'],
		#			'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' },
		#'Encelia_frutescens': {	'div': 0.03,
		#			'ind': ['CD102'],
		#			'seq': '/home/ssinghal/encelia/annotation/palmeri.annotated.fasta' }
	}	
			
np = 4
bowtie = '/home/ssinghal/bin/bowtie2'
bowtie_build = '/home/ssinghal/bin/bowtie2-build'
samtools = '/home/ssinghal/bin/samtools/samtools'
map_dir = '/home/ssinghal/encelia/variants/bam_files/'


def run_bowtie(map_dir, ind, seq, read_1, read_2, read_u, bowtie, np, samtools):
	o = open('run_bowtie.%s.sh' % ind, 'w')
	final = '%s%s.bowtie.palmeri_ref.bam' % (map_dir, ind)
	
	out1 = '%s%s1' % (map_dir, ind)
	out2 = '%s%s2' % (map_dir, ind)
	tmp = '%s%s.tmp.bam' % (map_dir, ind)
	finalout = final.replace('.bam', '')
		
	if not os.path.isfile(final):
		o.write("%s --mp 4 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p %s -x %s -1 %s -2 %s -S %s.sam -X 600\n" % (bowtie, np, seq, read_1, read_2, out1))
		o.write("%s --mp 4 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p %s -x %s -U %s -S %s.sam\n" % (bowtie, np, seq, read_u, out2))
		o.write("%s view -bS %s.sam > %s.bam\n" % (samtools, out1, out1))
		o.write("%s view -bS %s.sam > %s.bam\n" % (samtools, out2, out2))
		o.write("%s merge %s %s.bam %s.bam\n" % (samtools, tmp, out1, out2))
		o.write("%s sort %s %s\n" % (samtools, tmp, finalout))
		o.write("rm %s.sam %s.sam %s.gz %s.bam %s.bam %s\n" % (out1, out2, out1, out1, out2, tmp))
	o.close()
	return final


def run_smalt(map_dir, ind, seq, read_1, read_2, read_u, np, samtools):
        o = open('run_smalt.%s.sh' % ind, 'w')
        final = '%s%s.smalt.palmeri_ref.bam' % (map_dir, ind)

        out1 = '%s%s_smalt1.sam' % (map_dir, ind)
        out2 = '%s%s_smalt2.sam' % (map_dir, ind)
        tmp = '%s%s_smalt.tmp.bam' % (map_dir, ind)
        finalout = final.replace('.bam', '')

        if not os.path.isfile(final):
                o.write("smalt map -n 4 -o %s /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele %s %s\n" % (out1, read_1, read_2))
                o.write("smalt map -n 4 -o %s /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele %s\n" % (out2, read_u))
                o.write("%s view -bS %s > %s.bam\n" % (samtools, out1, out1))
                o.write("%s view -bS %s > %s.bam\n" % (samtools, out2, out2))
                o.write("%s merge %s %s.bam %s.bam\n" % (samtools, tmp, out1, out2))
                o.write("%s sort %s %s\n" % (samtools, tmp, finalout))
                o.write("rm %s %s %s.gz %s.bam %s.bam %s\n" % (out1, out2, out1, out1, out2, tmp))
        o.close()
        return final


def run_ngm(map_dir, ind, seq, read_1, read_2, read_u, np, samtools):
        o = open('run_ngm.%s.sh' % ind, 'w')
        final = '%s%s.ngm.palmeri_ref.bam' % (map_dir, ind)

        out1 = '%s%s_ngm1.sam' % (map_dir, ind)
        out2 = '%s%s_ngm2.sam' % (map_dir, ind)
        tmp = '%s%s_ngm.tmp.bam' % (map_dir, ind)
        finalout = final.replace('.bam', '')

        if not os.path.isfile(final):
                o.write("ngm -r %s -1 %s -2 %s -o %s -t 4\n" % (seq, read_1, read_2, out1))
                o.write("ngm -r %s -q %s -o %s -t 4\n" % (seq, read_u, out2))
                o.write("%s view -bS %s > %s.bam\n" % (samtools, out1, out1))
                o.write("%s view -bS %s > %s.bam\n" % (samtools, out2, out2))
                o.write("%s merge %s %s.bam %s.bam\n" % (samtools, tmp, out1, out2))
                o.write("%s sort %s %s\n" % (samtools, tmp, finalout))
                o.write("rm %s %s %s.gz %s.bam %s.bam %s\n" % (out1, out2, out1, out1, out2, tmp))
        o.close()
        return final


def run_samtools(seq, bams, samtools, bcftools_dir, map_dir, species):
	out_bcf = '%s%s.raw.bcf' % (map_dir, species)
	out_vcf = '%s%s.bowtie.samtools.vcf' % (map_dir, species)
	if not os.path.isfile(out_vcf):
		subprocess.call('%s mpileup -C 50 -uf %s %s | %sbcftools call -m -O b -v - > %s' % (samtools, seq, ' '.join(bams), bcftools_dir, out_bcf), shell=True)
		subprocess.call('%sbcftools view %s | %svcfutils.pl varFilter -w 0 > %s' % (bcftools_dir, out_bcf, bcftools_dir, out_vcf), shell=True)
	return out_vcf


for species in specieses:
	bams = []
	for ind in specieses[species]['ind']:
		seq = specieses[species]['seq']
		read_1 = '/home/ssinghal/encelia/trimmed_reads/%s/%s/%s_1_final.fastq' % (species, ind, ind)
		read_2 = read_1.replace('_1_', '_2_')
		read_u = read_1.replace('_1_', '_u_')

		# bam = run_bowtie(map_dir, ind, seq, read_1, read_2, read_u, bowtie, np, samtools)
		bam = run_smalt(map_dir, ind, seq, read_1, read_2, read_u, np, samtools)
		# bam = run_ngm(map_dir, ind, seq, read_1, read_2, read_u, np, samtools)
		#bams.append(bam)
	# species_vcf = run_samtools(seq, bams, samtools, bcftools_dir, map_dir, species)

