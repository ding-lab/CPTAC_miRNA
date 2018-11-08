###########################################################################################
### Annotation : Python Script
###########################################################################################
### @author: Sunantha Sethuraman (November, 2018)
###
### Description: Annotate miRNA sequencing reads using miRBase and GENCODE annotation information.
### The annotated information was downloaded and pre-processed and the script is not compatible with user
### supplied annotation, unless the information is pre-processed the same way.
###
### Purpose: The primary purpose and scope of this script is to processes CPTAC miRNA sequencing data.
###
### Usage: mirna_annotation.py -s <sample_name>
###########################################################################################

###########################################################################################
### Import modules
###########################################################################################

import pandas as pd
import re
import numpy as np
import string
import random
import argparse

###########################################################################################
### Parse arguments
###########################################################################################

parser = argparse.ArgumentParser(description='miRNA annotation pipeline')
parser.add_argument('--sample', metavar='s', type=str, nargs='+',
                   help='sample name')
                   
args = parser.parse_args()
s_name = args[0]
###########################################################################################
### Assign filenames and check if they exist
###########################################################################################

mirfile = s_name + "_mirbase.bed"
trnafile = s_name + "_tRNAs.bed"
exonfile = s_name + "_exons.bed"
unanfile = s_name + "_unannot.bed"
samfile = s_name + ".sam"

###########################################################################################
### Prioritization within miRNA annotations based on overlap
###########################################################################################

### All mature miRNA sequences will also be present in pre-miRNAs. Thus we need an extra over-lap based prioritization to distinguish mature miRNAs from pre-miRNAs. 
### Mature form is prioritized over pre form.

with open(mirfile, "r") as file:
	mirs = pd.read_csv(file, sep="\t", header=None)

mirs = mirs.iloc[:, [0,1,2,3,5,6,7,8,9,13,15,16]]
mirs.columns = ["Chr_read", "Start_read", "End_read", "Read_ID", "Strand", "Chr_ann", "Start_ann", "End_ann", "Annotation", "Class", "Tag", "Overlap" ]
mirs = mirs.sort_values(['Overlap', 'Class'], ascending = [False, True]).groupby(['Chr_read', 'Start_read', 'End_read', 'Read_ID']).head(1)

###########################################################################################
### tRNA annotation clean up
###########################################################################################

with open(trnafile, "r") as file:
	trna = pd.read_csv(file, sep = "\t", header = None)

trna = trna.iloc[:, [0,1,2,3,5,6,7,8,9,13,15,16]]
trna.columns = ["Chr_read", "Start_read", "End_read", "Read_ID", "Strand", "Chr_ann", "Start_ann", "End_ann", "Annotation", "Class", "Tag", "Overlap" ]
trna.Annotation = [re.sub(".*gene_type=(.*?);.*", r"\1", x) for x in trna.Tag]

###########################################################################################
### Exon annotation clean up
###########################################################################################

with open(exonfile, "r") as file:
	exon = pd.read_csv(file, sep = "\t", header = None)

exon = exon.iloc[:, [0,1,2,3,5,6,7,8,9,13,15,16]]
exon.columns = ["Chr_read", "Start_read", "End_read", "Read_ID", "Strand", "Chr_ann", "Start_ann", "End_ann", "Annotation", "Class", "Tag", "Overlap" ]
exon.Class = [re.sub(".*gene_type=(.*?);.*", r"\1", x) for x in exon.Tag]

###########################################################################################
### Unannotated reads clean up
###########################################################################################

with open(unanfile, "r") as file:
	unan = pd.read_csv(file, sep="\t", header=None)


unan = unan.iloc[:, [0,1,2,3,5]]
unan.columns = ["Chr_read", "Start_read", "End_read", "Read_ID", "Strand"]
unan = unan.reindex(columns=["Chr_read", "Start_read", "End_read", "Read_ID", "Strand", "Chr_ann", "Start_ann", "End_ann", "Annotation", "Class", "Tag", "Overlap" ], fill_value=".")
unan.Class = "unannotated"

###########################################################################################
### Combine all annotations
###########################################################################################

comb = mirs.append(trna)
comb = comb.append(exon)
comb = comb.append(unan)

###########################################################################################
### Remove multiple annotations of the same mapping based on a set priority
###########################################################################################

rna_priorities = ['miRNA', 'miRNA_primary_transcript', 'tRNA', 'rRNA', 'sRNA', 'snRNA', 'snoRNA', 'scaRNA', 'scRNA', 'ribozyme', 'vaultRNA', 'Mt_tRNA',  'Mt_rRNA', 'protein_coding', 
					'3prime_overlapping_ncRNA', 'bidirectional_promoter_lncRNA', 'lincRNA','macro_lncRNA', 'misc_RNA', 'non_coding',
					'antisense', 'sense_intronic', 'sense_overlapping',
					'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene', 'IG_C_pseudogene','IG_J_pseudogene', 'IG_V_pseudogene', 'IG_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene',
					'rRNA_pseudogene', 'processed_transcript', 'pseudogene', 'polymorphic_pseudogene', 'processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 
					'unitary_pseudogene', 'unprocessed_pseudogene', 'TEC','unannotated']
					
comb['Class'] = pd.Categorical(comb['Class'], rna_priorities)
comb = comb.sort_values('Class', ascending = True).groupby(['Chr_read', 'Start_read', 'End_read', 'Read_ID']).head(1)

ann = comb.copy() # Save the annotation information at this stage for later use in annotating the SAM file

###########################################################################################
### Multi-mapping hit management : Multi-mappers within same annotation class
###########################################################################################

## For multi-mapping hits with the same annotation class, pick one random hit to keep
## The code below is in some sense equivalent to running : comb = comb.groupby(['Read_ID', 'Class']).apply(lambda x: x.sample(n=1))
## However this is computationally intensive, so we use the lengthy code below. 

randomNum = range(0, comb.shape[0])
random.shuffle(randomNum)
comb['randomNum'] = randomNum
mhits = comb[comb['Read_ID'].groupby(comb['Read_ID']).transform('size') >1]
mhits = mhits.sort_values('randomNum', ascending = True).groupby(['Read_ID', 'Class']).head(1)
comb = comb.loc[~ comb['Read_ID'].isin(mhits.Read_ID)] 
comb = comb.append(mhits, ignore_index = True)
comb = comb.drop('randomNum', axis =1)

###########################################################################################
### Multi-mapping hit management : Multi-mappers across annotation classes
###########################################################################################

## For multi-mapping hits across annotation classes, use the RNA priorities as defined before

comb['Class'] = pd.Categorical(comb['Class'], rna_priorities) 	# Not necessary, better to ensure
comb = comb.sort_values('Class', ascending = True).groupby('Read_ID').head(1)

###########################################################################################
### Write annotation output as a TXT file
###########################################################################################
### In this output, each mapped read, has a single annotation. Multiple mapping and multiple annotations are not reflected in this file.
### For multiple-mapping annotation, see annotated SAM file generated below where are mapping is preserved.

comb.to_csv(s_name + "_Annotated.txt", sep='\t', index=False)

###########################################################################################
### SAM file annotation
###########################################################################################

### Format annotation information for SAM annotation

ann = ann[['Read_ID', 'Chr_read', 'Start_read', 'Annotation', 'Class', 'Tag']]
ann.Class = ann.Class.astype(str)
ann.Start_read = ann.Start_read + 1
ann.Annotation = "YA:Z:"+ ann.Annotation 	# This structure is essential to make this information readable by samtools
ann.Class = "YC:Z:"+ ann.Class
ann.Tag = "YT:Z:" + ann.Tag
ann = ann.astype('str')	  # Pandas merge runs into problems unless dtypes match, this ensures that

### Read SAM file and format for annotation

with open(samfile, "r") as file:
	sam = pd.read_csv(file, sep = "\t", names = list(string.ascii_lowercase)[0:20], engine = 'c', comment='@', dtype = 'str')

sam = sam.rename(index=str, columns={"a": "Read_ID", "c": "Chr_read", "d": "Start_read"})
sam = sam.astype('str')

### Annotation

sam_ann = pd.merge(sam, ann, 'left')
sam_ann = sam_ann.replace("nan", '', regex=True)

cols = sam_ann.columns.tolist()		# SAM format does not support empty columns, this this rearrangement is essential
cols = cols[0:11] + cols[-3:] + cols[11:-3]
sam_ann = sam_ann[cols]

###########################################################################################
### Write annotation output as a SAM file
###########################################################################################
### In this output, each line represents a unique mapping. Some reads might have multiple (up to 10) lines (mappings).
### The output is NOT ready for viewing using samtools. The pipeline will manipulate this further to make it samtools compatible.

sam_ann.to_csv(s_name + "_Annotated.sam", sep='\t', index=False, header=False)
