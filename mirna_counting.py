###########################################################################################
### Annotation : Python Script
###########################################################################################
### @author: Sunantha Sethuraman (November, 2018)
###
### Description: Count miRNA reads using the output of mirna_annotation.py. This script will write three output files
### 1. <sample_name>_mature.txt: Counts of mature isoforms of miRNAs
### 2. <sample_name>_precursor.txt: Counts of precursor (stem-loop) form of miRNAs
### 3. <sample_name>_miRNA_expression.txt: Combined counts of mature and precursor forms of miRNAs.
###
### Purpose: The primary purpose and scope of this script is to processes CPTAC miRNA sequencing data.
###
### Usage: mirna_counting.py --sample <sample_name>  --mirbase <mirbase_file>
###########################################################################################

###########################################################################################
### Import modules
###########################################################################################

import pandas as pd
import re
import argparse

###########################################################################################
### Parse arguments
###########################################################################################

parser = argparse.ArgumentParser(description='miRNA annotation pipeline')
parser.add_argument('--sample', type=str, help='sample name', required = True)
parser.add_argument('--mirbase', type=str, help='mirbase file', required = True)
                   
args = parser.parse_args()
s_name = args.sample
m_file = args.mirbase

###########################################################################################
### Read input files
###########################################################################################

annfile = s_name + "_Annotated.txt"

with open(annfile, "r") as file:
	ann = pd.read_csv(file, sep="\t")

with open(m_file, "r") as file:
	mbase = pd.read_csv(file, sep="\t", header= None)

###########################################################################################
### Clean up miRBase file (to be used later)
###########################################################################################

mbase = mbase.loc[mbase[7] == "miRNA_primary_transcript", :]
tag = pd.Series(re.sub(".*?=(.*;).*?=(.*;).*?=(.*)", r"\1\2\3", x) for x in mbase[9])
mbase = pd.DataFrame(tag.str.split(';').tolist(), columns = ['ID','Alias', 'Name'])

###########################################################################################
### Count mature miRNAs
###########################################################################################

isomirs = ann.loc[ann['Class'] == "miRNA", :]

### Ensure only miRBase annotations are present (Not necessary)
cond = isomirs['Tag'].str.contains('ID=MIMAT')
isomirs = isomirs.loc[cond,:]

### Save the number of reads aligned to mature miRNAs
mir_reads = isomirs.shape[0]

### Count
tag = pd.Series(re.sub(".*?=(.*;).*?=(.*;).*?=(.*;).*?=(.*)", r"\1\2\3\4", x) for x in isomirs.Tag)
reads = pd.DataFrame(tag.str.split(';').tolist(), columns = ['ID','Alias', 'Name', 'Derives_from'])
reads = reads.drop_duplicates()
isomirs['miRNA_Name'] = [re.sub(".*Name=(.*?);.*", r"\1", x) for x in isomirs.Tag]
counts = isomirs[['miRNA_Name', 'Annotation']].groupby('Annotation').count()
counts.columns = ['miR_Count_Raw']
mir_counts = pd.merge(reads, counts, left_on = 'ID', right_index=True)

###########################################################################################
### Count precursor miRNAs
###########################################################################################

premirs = ann.loc[ann['Class'] == "miRNA_primary_transcript", :]

### Ensure only miRBase annotations are present (Not necessary)
cond = premirs['Tag'].str.contains('ID=MI0')
premirs = premirs.loc[cond,:]

### Save the number of reads aligned to precursor miRNAs
pre_reads = premirs.shape[0]

### Count
tag = pd.Series(re.sub(".*?=(.*;).*?=(.*;).*?=(.*)", r"\1\2\3", x) for x in premirs.Tag)
reads = pd.DataFrame(tag.str.split(';').tolist(), columns = ['ID','Alias', 'Name'])
reads = reads.drop_duplicates()
premirs['primary_miRNA_Name'] = [re.sub(".*Name=(.*?);.*", r"\1", x) for x in premirs.Tag]
counts = premirs[['primary_miRNA_Name', 'Annotation']].groupby('Annotation').count()
counts.columns = ['primary_miR_Count_Raw']
pre_counts = pd.merge(reads, counts, left_on = 'ID', right_index=True)

###########################################################################################
### Count total miRNAs
###########################################################################################

mats = pd.DataFrame(mir_counts[['Derives_from', 'miR_Count_Raw']])
pres = pd.DataFrame(pre_counts[['ID', 'primary_miR_Count_Raw']])
mats.columns = ['ID', 'Count_Raw']
pres.columns = ['ID', 'Count_Raw']

all = mats.append(pres, ignore_index = True)
counts = all.groupby('ID').sum()
all_counts = pd.merge(mbase, counts, left_on = 'ID', right_index=True, how = 'right')

###########################################################################################
### Calculate TPM
###########################################################################################

### Transcripts per million reads (TPM) is calculated as: (number of reads of each miRNA * 1000000) / (total reads)
### Here "total reads" represents the total number of reads aligned to either a mature or precursor miRNA sequence
### It is important to note that the "total reads" is not the same the library size / sequencing depth.

total_reads = mir_reads + pre_reads
mir_counts['TPM'] = mir_counts['miR_Count_Raw']*1000000 / total_reads
pre_counts['TPM'] = pre_counts['primary_miR_Count_Raw']*1000000 / total_reads
all_counts['TPM'] = all_counts['Count_Raw']*1000000/ total_reads

### Write output files 

mir_counts.to_csv(s_name + "_mature.txt", sep='\t', index=False)
pre_counts.to_csv(s_name + "_precursor.txt", sep='\t', index=False)
all_counts.to_csv(s_name + "_miRNA_expression.txt", sep = "\t", index = False)
