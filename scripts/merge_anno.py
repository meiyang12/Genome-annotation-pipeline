from sys import argv

seq_list=[x.strip().split()[0][1:] for x in open(argv[1])]
uniportGo={x.split('\t')[0]:x.strip().split('\t')[1] for x in open(argv[2])}
kegg={x.split('\t')[0]:x.strip().split('\t')[1] for x in open(argv[3])}
pfam={x.split('\t')[0]:x.strip().split('\t')[1] for x in open(argv[4])}
for x in seq_list:
    uni=''
    ke=''
    pf=''
    if x in uniportGo:
        uni=uniportGo[x]
    if x in kegg:
        ke=kegg[x]
    if x in pfam:
        pf=pfam[x]
    print(x+'\t'+uni+'\t'+ke+'\t'+pf)