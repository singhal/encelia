import re
import subprocess
import glob
import random

specieses = ['asperifolia', 'canescens', 'farinosa', 'frutescens', 'palmeri', 'ventorum', 'helianthus']
genes = {}
out_dir =  '/home/ssinghal/encelia/analysis/gene_trees/'

def get_seqs(species, genes):
	f = open('/home/ssinghal/encelia/annotation/%s.annotated.fasta' % species, 'r')
	for l in f:
		if re.search('>', l):
			d = re.split('\t', l.rstrip())
			if len(d) > 1:
				gene = d[2]
				seq = f.next().rstrip()
				cds_start = int(re.search('gs(\d+)', d[1]).group(1))
				cds_end = int(re.search('ge(\d+)', d[1]).group(1))

				cds = seq[(cds_start - 1):cds_end]

				if gene not in genes:
					genes[gene] = {}
				genes[gene][species] = cds
	return genes


def write_files(out_dir, genes):
	to_align_files = []
	for gene in genes:
		# all the species are here!
		if len(genes[gene].keys()) == len(specieses):
			out_file = '%s/alignments/%s.fa' % (out_dir, gene)
			o = open(out_file, 'w')
			for species in genes[gene]:
				o.write('>%s\n%s\n' % (species, genes[gene][species]))
			o.close()
			to_align_files.append(out_file)
	return to_align_files	


def align_files(to_align_files):
	aligned_files = []
	for file in to_align_files:
		aln = file + '.aln'
		subprocess.call('muscle -in %s -out %s' % (file, aln), shell=True)
		aligned_files.append(aln)
	return aligned_files


def convert_phyml(aligned_files):
	new = []
	for file in aligned_files:
		out = file.replace('.fa.aln', '.aln.phy')
		f = open(file, 'r')
		seq = {}
		id = ''
		for l in f:
			if re.search('>', l):
				id = l.rstrip()
				id = id.replace('>', '')
				seq[id] = ''
			else:
				seq[id] += l.rstrip()
		f.close()

		o = open(out, 'w')
		o.write(' %s %s\n' % (len(seq), len(seq.values()[0])))
		for id in seq:
			id2 = id + (20 - len(id)) * ' '
			o.write('%s%s\n' % (id2, seq[id]))
		o.close()
		new.append(out)
	return new


def make_trees(out_dir, aligned_files):
	out_dir = out_dir + 'trees/'
	for file in aligned_files:
		name = re.search('(Migut\S+).aln', file).group(1)
		rand_num = random.randint(1,10000)
		rand_num2 = random.randint(1,10000)
		subprocess.call('raxmlHPC -f a -o helianthus -m GTRCAT -n %s -s %s -w %s -N 100 -x %s -p %s' % (name, file, out_dir, rand_num, rand_num2), shell=True)
	

#for species in specieses:
#	genes = get_seqs(species, genes)
#to_align_files = write_files(out_dir, genes)
#aligned_files = align_files(to_align_files)
aligned_files = glob.glob('/home/ssinghal/encelia/analysis/gene_trees/alignments/*phy')
#aligned_files = convert_phyml(aligned_files)
make_trees(out_dir, aligned_files)	
