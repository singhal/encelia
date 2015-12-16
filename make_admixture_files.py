import re
import gzip
import random

vcf1 = '/home/ssinghal/encelia/variants/vcfs/palmeri.ngm.annotated.filtered_qual.high_cov.filtered_fs.vcf.gz'
vcf2 = '/home/ssinghal/encelia/variants/vcfs/ventorum.ngm.annotated.filtered_qual.high_cov.filtered_fs.vcf.gz'
ind = {'palmeri': 5, 'ventorum': 5}

snps = {}

def get_vcf(vcf, snps, sp):
	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l)
			if not re.search(',', d[4]):
				if d[0] not in snps:
					snps[d[0]] = {}
				if d[1] not in snps[d[0]]:
					snps[d[0]][d[1]] = {'palmeri': [], 'ventorum': [], 'alt': d[4]}
				if d[4] == snps[d[0]][d[1]]['alt']:
					for geno in d[9:]:
						if re.match('0/0', geno):
							snps[d[0]][d[1]][sp].append(0)
						elif re.match('0/1', geno):
							snps[d[0]][d[1]][sp].append(1)
						elif re.match('1/1', geno):
							snps[d[0]][d[1]][sp].append(2)
						else:
							snps[d[0]][d[1]][sp].append(9)
	f.close()
	return snps

snps = get_vcf(vcf1, snps, 'palmeri')
snps = get_vcf(vcf2, snps, 'ventorum')


for chr in snps:
        for pos in snps[chr]:
                for sp in ind:
                        if len(snps[chr][pos][sp]) == 0:
                                snps[chr][pos][sp] = [0] * ind[sp]


def get_blocks(values, dist):
    mi, ma = 0, 0
    result = []
    temp = []
    for v in sorted(values):
        if not temp:
            mi = ma = v
            temp.append(v)
        else:
            if abs(v - mi) < dist and abs(v - ma) < dist:
                temp.append(v)
                if v < mi:
                    mi = v
                elif v > ma:
                    ma = v
            else:
                if len(temp) > 0:
                    result.append(temp)
                mi = ma = v
                temp = [v]
    return result


out = open('/home/ssinghal/encelia/analysis/admixture/encelia.admixture.geno', 'w')
for chr in snps:
	positions = sorted([int(x) for x in snps[chr].keys()])
	positions = get_blocks(positions, 100)
	positions = [str(random.choice(x)) for x in positions]	

        for pos in positions:
                snps_str = ''
                for species in ['palmeri', 'ventorum']:
                        for bp in snps[chr][pos][species]:
                                snps_str += str(bp)
                out.write(snps_str + '\n')
	if len(positions) > 0:
		print chr
out.close()
