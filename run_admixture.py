import subprocess

for k in range(1,11):
	print '/home/ssinghal/bin/admixture_linux-1.23/admixture --cv /home/ssinghal/encelia/analysis/admixture/encelia.admixture.geno %s > log%s.out' % (k, k)
