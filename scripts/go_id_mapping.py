from sys import argv

swissprot_out = {y.split()[0]: y.split()[0].split('|')[1]
                 for y in open(argv[2])}

id_hash = {}
for x in open(argv[1]):
    tmp = x.split('\t')
    tmp_key = tmp[0]
    tmp_value = str('')
    for y in tmp:
        if 'GO' in y:
            tmp_value = tmp_value+str(y)
    id_hash[tmp_key] = tmp_value

seq_dict = {}
for x in open(argv[3]):
    if x.startswith('>'):
        seq_dict[x.split()[0][1:]] = x.split()[1].split(':')[1].strip()

for x in seq_dict:
    if x in swissprot_out:	
        if swissprot_out[x] in id_hash:
            print(x+'\t'+seq_dict[x]+'\t'+swissprot_out[x]+'\t'+id_hash[swissprot_out[x]].replace(';', ','))
        else:
            print(x+'\t'+seq_dict[x]+'\t'+swissprot_out[x]+'\t')
    else:
        print(x+'\t'+seq_dict[x]+'\t')
    
