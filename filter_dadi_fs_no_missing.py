import re

original = '/home/ssinghal/encelia/analysis/dadi/palmeri_ventorum.same_ref.fs'
new = '/home/ssinghal/encelia/analysis/dadi/palmeri_ventorum.same_ref.no_missing.fs'

f = open(original, 'r')
o = open(new, 'w')

head = f.next().rstrip()
o.write(head + '\n')

for l in f:
	d = re.split('\s+', l.rstrip())
	if (int(d[3]) + int(d[6])) == 10 and (int(d[4]) + int(d[7])) == 10:
		o.write(l)
f.close()
o.close()
