import glob
import subprocess
import os

files = glob.glob("/home/ssinghal/encelia/variants/bam_files/*bam")
for file in files:
	if not os.path.isfile(file + '.bai'):
		subprocess.call("samtools index %s" % (file), shell=True)

