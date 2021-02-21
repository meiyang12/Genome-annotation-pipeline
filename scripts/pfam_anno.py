from sys import argv

hmm_dict = {}
for x in open(argv[1]):
    if not x.startswith('#'):
        tmp = x.strip().split(maxsplit=18)
        if tmp[2] not in hmm_dict:
            hmm_dict[tmp[2]] = []
            hmm_dict[tmp[2]].append([tmp[1], tmp[0], tmp[-1]])
        else:
            hmm_dict[tmp[2]].append([tmp[1], tmp[0], tmp[-1]])

for x in hmm_dict:
    pf = ''
    for y in hmm_dict[x]:
        pf += y[0]+'('+y[1]+'('+y[2]+'));'
    print(x+'\t'+pf)
