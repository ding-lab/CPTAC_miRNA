###########################################################################################
### @Author: Sunantha Sethuraman (s.sethuraman@wustl.edu)
### @Date: November 2018
### @Description: Annotation information download and clean up
### @Availability: Annotation files as used in the CPTAC processing pipeline are provided here on github
###########################################################################################


###########################################################################################
### Downloaded GENCODE v29 BASIC annotation
### Active link at the time of download: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gff3.gz
###########################################################################################

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gff3.gz
 
gunzip gencode.v29.basic.annotation.gff3.gz

gff2bed < gencode.v29.chr_patch_hapl_scaff.basic.annotation.gff3 >  gencode.v29.chr_patch_hapl_scaff.basic.annotation.bed

grep "^chr" gencode.v29.chr_patch_hapl_scaff.basic.annotation.bed > gencode_std.bed

grep "ID=exon" gencode_std.bed | grep -v "gene_type=miRNA" > gencode_std_exons.bed


###########################################################################################
### Downloaded miRBase v22 annotation
### Active link at the time of download: ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
###########################################################################################

wget -c ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3

gff2bed < hsa.gff3 > hsa.bed

cp hsa.bed mirbase.bed


###########################################################################################
### Downloaded GENCODE v29 tRNA annotation
### Active link at the time of download: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.tRNAs.gff3.gz
###########################################################################################

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.tRNAs.gff3.gz

gunzip gencode.v29.tRNAs.gff3.gz

gff2bed < gencode.v29.tRNAs.gff3 > gencode.v29.tRNAs.bed

cp gencode.v29.tRNAs.bed gencode_tRNAs.bed

###########################################################################################
