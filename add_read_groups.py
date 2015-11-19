import re
import subprocess
import glob
import os

bam_files = glob.glob('/home/ssinghal/encelia/variants/bam_files/*bam')
java = '/home/ssinghal/bin/jre1.7.0_67/bin/java'
picard_dir = '/home/ssinghal/bin/picard-tools/'

for bam_file in bam_files:
	out = bam_file.replace('.bam', '.rg.bam')
	lib = re.search('(CD\d+)', bam_file).group(1)
	if not os.path.isfile(out):
		subprocess.call("%s -Xmx6g -jar %sAddOrReplaceReadGroups.jar I=%s O=%s SO=coordinate RGID=Encelia RGLB=%s RGPL=Illumina RGPU=GAIIx RGSM=%s MAX_RECORDS_IN_RAM=1000000 TMP_DIR=/home/ssinghal/tmp/" % (java, picard_dir, bam_file, out, lib, lib), shell=True)
		subprocess.call("samtools index %s" % (out), shell=True)
