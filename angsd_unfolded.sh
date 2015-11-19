# identify the sites that are variable in either lineage
# -P is the number of threads
/home/ssinghal/bin/angsd0.612/angsd -GL 1 -b /home/ssinghal/encelia/analysis/angsd/ventorum.bamlist.txt -anc /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa -rf /home/ssinghal/encelia/analysis/angsd/coverage.angsd_regions.txt -P 12 -out /home/ssinghal/encelia/analysis/angsd/ventorum -doSaf 1 -only_proper_pairs 0
/home/ssinghal/bin/angsd0.612/angsd -GL 1 -b /home/ssinghal/encelia/analysis/angsd/palmeri.bamlist.txt -anc /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa -rf /home/ssinghal/encelia/analysis/angsd/coverage.angsd_regions.txt -P 12 -out /home/ssinghal/encelia/analysis/angsd/palmeri -doSaf 1 -only_proper_pairs 0
# get the SFS for each pop
# first argument is saf file, second argument is the number of chromosomes, -P 24 is the number of cores we want to use
/home/ssinghal/bin/angsd0.612/realSFS /home/ssinghal/encelia/analysis/angsd/ventorum.saf 10 -P 12 > /home/ssinghal/encelia/analysis/angsd/ventorum.saf.sfs
/home/ssinghal/bin/angsd0.612/realSFS /home/ssinghal/encelia/analysis/angsd/palmeri.saf 10 -P 12 > /home/ssinghal/encelia/analysis/angsd/palmeri.saf.sfs
# want to identify sites that are variable in either of the populations
gunzip -c /home/ssinghal/encelia/analysis/angsd/ventorum.saf.pos.gz /home/ssinghal/encelia/analysis/angsd/palmeri.saf.pos.gz |sort  -S 50\%|uniq -d|sort -k1,1  -S 50\% > /home/ssinghal/encelia/analysis/angsd/palmeri_ventorum.intersect_sites.txt
# want to estimate genotypes at variable sites for both chromosomes
/home/ssinghal/bin/angsd0.612/angsd -GL 1 -b /home/ssinghal/encelia/analysis/angsd/ventorum.bamlist.txt -anc /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa -rf /home/ssinghal/encelia/analysis/angsd/coverage.angsd_regions.txt -P 12 -out /home/ssinghal/encelia/analysis/angsd/ventorum2 -sites /home/ssinghal/encelia/analysis/angsd/palmeri_ventorum.intersect_sites.txt -doSaf 1 -only_proper_pairs 0
/home/ssinghal/bin/angsd0.612/angsd -GL 1 -b /home/ssinghal/encelia/analysis/angsd/palmeri.bamlist.txt -anc /home/ssinghal/encelia/variants/ancestral_allele/ancestral_allele.fa -rf /home/ssinghal/encelia/analysis/angsd/coverage.angsd_regions.txt -P 12 -out /home/ssinghal/encelia/analysis/angsd/palmeri2 -sites /home/ssinghal/encelia/analysis/angsd/palmeri_ventorum.intersect_sites.txt -doSaf 1 -only_proper_pairs 0
# estimate joint frequency
/home/ssinghal/bin/angsd0.612/realSFS 2dsfs /home/ssinghal/encelia/analysis/angsd/ventorum2.saf /home/ssinghal/encelia/analysis/angsd/palmeri2.saf 10 10 -P 12 > /home/ssinghal/encelia/analysis/ventorum_palmeri.2D.sfs
