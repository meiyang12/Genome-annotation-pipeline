from sys import argv

swiss_dict = {}
for x in open(argv[1]):
    if x.startswith('>'):
        tmp = x.strip().split(' ', 1)
        key = tmp[0].split('|')[1]
        des = tmp[1]
        swiss_dict[key] = des

seq_dict = {}
for x in open(argv[2]):
    x = x.strip()
    if x.startswith('>'):
        key = x[1:]
        seq_dict[key] = []
    else:
        seq_dict[key].append(x)

blast_dict = {}
for x in open(argv[3]):
    tmp = x.split()
    seq_id = tmp[0]
    swiss_id = tmp[1].split('|')[1]
    score = float(tmp[-1])
    if seq_id not in blast_dict:
        blast_dict[seq_id] = [swiss_id, score]
    else:
        if score >= blast_dict[seq_id][1]:
            blast_dict[seq_id] = [swiss_id, score]
        else:
            pass

for x in seq_dict:
    if x not in blast_dict:
        print('>'+x+' annotation: Putative uncharacterized protein\n' +
              ''.join(seq_dict[x]))
    else:
        print('>'+x+' annotation: ' +
              swiss_dict[blast_dict[x][0]]+'\n'+''.join(seq_dict[x]))
