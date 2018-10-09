import itertools

keywords = [''.join(i) for i in itertools.product(['A','C','G','T'], repeat = 4)]
new_keywords = [">Adapter_" + str(i) + "\n" + keywords[i]+ "TGGAATTCTCGGGTGCCAAGGC" for i in range(len(keywords))]

with open("/diskmnt/Projects/Users/ssethura/cptac_mirna/FASTQ/Adapters.fa", 'w') as f:
    f.write("\n".join(map(str, new_keywords)))