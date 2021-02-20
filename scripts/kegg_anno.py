from sys import argv

kegg_dict={x.split()[0][1:]:x.split()[1] for x in open(argv[1]) if x.startswith('>')}
diamond_dict={x.split()[0]:x.split()[1] for x in open(argv[2])}
for x in diamond_dict:
	print(x+'\t'+kegg_dict[diamond_dict[x]])
