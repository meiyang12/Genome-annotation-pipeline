from sys import argv
import re

with open(argv[1]) as f:
    for x in f:
        x=x.strip()
        if not x.startswith('#'):
            if 'mRNA' in x:
                n=1
                print(x)
            if 'CDS' in x:
                tmp=x.split('\t')
                id=re.findall(r'Parent=(.*?);', tmp[-1])[0]+'.'+str(n)
                n+=1
                tmp[-1]='ID='+id+';'+tmp[-1]
                print('\t'.join(tmp))
