import re
import subprocess
import glob
import random
import numpy as np
import os

specieses = ['asperifolia', 'canescens', 'farinosa', 'frutescens', 'palmeri', 'ventorum']
genes = {}
out_dir =  '/home/ssinghal/encelia/analysis/dnds/'
prot_db =  '/home/ssinghal/encelia/genomes/Mguttatus_v2.0_256.fasta'


def get_protein(prot_db):
	protein = {}
	id = ''
	f = open(prot_db, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			protein[id] = ''
		else:
			seq = l.rstrip()
			seq = seq.replace('*', '')
			protein[id] += seq
	f.close()
	return protein


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


def translate(seq):
	map = {	'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    		'TCT':'S', 'TCC':'s', 'TCA':'S', 'TCG':'S',
  		'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    		'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    		'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    		'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    		'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    		'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
   		'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    		'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    		'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    		'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    		'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    		'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    		'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

	triplets = [seq[i:i+3] for i in range(0, len(seq), 3)]
	aa = ''
	for triplet in triplets:
		if triplet in map:
			aa += map[triplet]
		else:
			aa += 'X'
	return aa


def check_frame(prots, seq, gene):
	for i in range(0,3):
		# translate that frame
		aa = translate(seq[i:])
		# do the blast
		query = 'query.fa'
		q = open(query, 'w')
		q.write('>query\n%s\n' % aa)
		q.close()
		subject = 'subject.fa'
                s = open(subject, 'w')
                s.write('>subject\n%s\n' % prots[gene])
                s.close()

		blast = subprocess.Popen('blastp -query %s -subject %s -outfmt 6' % (query, subject), shell=True, stdout=subprocess.PIPE)
		top_match = 0 
		for l in blast.stdout:
			d = re.split('\t', l.rstrip())
			if float(d[2]) > top_match:
				top_match = float(d[2])
		
		if (top_match > 30):
			return i
	return None
			

def write_files(prots, out_dir, genes):
	to_align_files = []
	for gene in genes:
		# are the focal species in this gene?
		if 'ventorum' in genes[gene] and 'palmeri' in genes[gene]:
			max_length = np.max([len(x) for x in genes[gene].values()])
			out_file = '%s/alignments/%s.fa' % (out_dir, gene)
			o = open(out_file, 'w')
			for species in genes[gene]:
				# check frame
				frame = check_frame(prots, genes[gene][species], gene)
				# check length
				if frame != None:
					print frame
					aa = translate(genes[gene][species][frame:])
					len_seq = len(re.search('^([^\*]+)', aa).group(1)) * 3
					if (len_seq / float(max_length)) > 0.8:
						seq = genes[gene][species][frame:(frame+len_seq)]
						o.write('>%s\n%s\n' % (species, seq))
			o.close()
			to_align_files.append(out_file)
	return to_align_files	


def align_files(to_align_files):
	aligned_files = []
	for file in to_align_files:
		prot_file = file.replace('.fa', '.prot.fa')
		o = open(prot_file, 'w')
		f = open(file, 'r')
		for l in f:
			if re.search('>', l):
				o.write(l)
			else:
				aa = translate(l.rstrip())
				o.write(aa + '\n')
		o.close()
		f.close()
		aln = prot_file + '.aln'
		subprocess.call('muscle -in %s -out %s' % (prot_file, aln), shell=True)
		aligned_files.append(aln)
	return aligned_files


def gapsFromPeptide( peptide_seq, nucleotide_seq ):
    """ Transfers gaps from aligned peptide seq into codon partitioned nucleotide seq (codon alignment) 
          - peptide_seq is an aligned peptide sequence with gaps that need to be transferred to nucleotide seq
          - nucleotide_seq is an un-aligned dna sequence whose codons translate to peptide seq"""
    def chunks(l, n):
        """ Yield successive n-sized chunks from l."""
        for i in xrange(0, len(l), n):
            yield l[i:i+n]
    codons = [codon for codon in chunks(nucleotide_seq,3)]  #splits nucleotides into codons (triplets) 
    gappedCodons = []
    codonCount = 0
    for aa in peptide_seq:  #adds '---' gaps to nucleotide seq corresponding to peptide
        if aa!='-':
            gappedCodons.append(codons[codonCount])
            codonCount += 1
        else:
            gappedCodons.append('---')
    return(''.join(gappedCodons))


def convert_phyml(aligned_files):
	new = []
	for aa_file in aligned_files:
		out = aa_file.replace('.prot.fa.aln', '.aln.phy')
				
		f = open(aa_file, 'r')
		aa = {}
		id = ''
		for l in f:
			if re.search('>', l):
				id = l.rstrip()
				id = id.replace('>', '')
				aa[id] = ''
			else:
				aa[id] += l.rstrip()
		f.close()
	
		nuc_file = aa_file.replace('.prot.fa.aln', '.fa')
		f = open(nuc_file, 'r')
                nuc = {}
                id = ''
                for l in f:
                        if re.search('>', l):
                                id = l.rstrip()
                                id = id.replace('>', '')
                                nuc[id] = ''
                        else:
                                nuc[id] += l.rstrip()
                f.close()


		if len(aa) > 1:
			o = open(out, 'w')
			o.write(' %s %s\n' % (len(aa), 3*len(aa.values()[0])))
			for id in aa:
				id2 = id + (20 - len(id)) * ' '
				nuc_seq = gapsFromPeptide(aa[id], nuc[id])
				o.write('%s%s\n' % (id2, nuc_seq))
			o.close()
			new.append(out)
	return new


def make_trees(out_dir, aligned_files):
        out_dir = out_dir + 'trees/'
        for file in aligned_files:
                name = re.search('(Migut\S+).aln', file).group(1)
                rand_num = random.randint(1,10000)
                rand_num2 = random.randint(1,10000)
                subprocess.call('raxmlHPC -f a -m GTRCAT -n %s -s %s -w %s -N 100 -x %s -p %s' % (name, file, out_dir, rand_num, rand_num2), shell=True)


def run_paml(out_dir, aligned_files):
	constant_seq = '%sseq.phy' % out_dir
	constant_tree = '%sseq.tre' % out_dir
	constant_out = '%sout' % out_dir

	for file in aligned_files:
		name = re.search('(Migut.*p)\.aln', file).group(1)
		tree = '%strees/RAxML_bestTree.%s' % (out_dir, name)
		out = '%sresults/%s.paml.out' % (out_dir, name)

		if os.path.isfile(tree):
			subprocess.call('cp %s %s' % (file, constant_seq), shell=True)
			subprocess.call('cp %s %s' % (tree, constant_tree), shell=True)
			subprocess.call('~/bin/paml4.8/bin/codeml ~/bin/paml4.8/bin/codeml.ctl', shell=True)
			subprocess.call('cp %s %s' % (constant_out, out), shell=True)

#prots = get_protein(prot_db)
#for species in specieses:
#	genes = get_seqs(species, genes)
# to_align_files = write_files(prots, out_dir, genes)
aligned_files = glob.glob('/home/ssinghal/encelia/analysis/dnds/alignments/*aln')
# aligned_files = align_files(to_align_files)
aligned_files = convert_phyml(aligned_files)
# make_trees(out_dir, aligned_files)
run_paml(out_dir, aligned_files)
