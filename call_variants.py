import re

species = {	'palmeri': ['CD659', 'CD661', 'CD662', 'CD663', 'CD664'],
		'ventorum': ['CD011', 'CD012', 'CD013', 'CD014', 'CD015']}

for sp in species:
	for type in ['bowtie', 'ngm', 'smalt']:
		o = open('%s_%s.sh' % (sp, type), 'w')
		bam_files = ['%s.%s.palmeri_ref.bam' % (ind, type) for ind in species[sp]]
		o.write('#~/bin/samtools/samtools mpileup -A -I -ugf /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa -o %s.%s.raw.bcf %s\n' % (sp, type, ' '.join(bam_files)))
		o.write('~/bin/bcftools/bcftools call -vmO z -o %s.%s.raw.vcf.gz %s.%s.raw.bcf\n' % (sp, type, sp, type))
		o.close()
