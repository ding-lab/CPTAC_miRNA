###########################################################################################
### Activate conda environment (ignore if all the necessary softwares are installed in the local environment)
###########################################################################################

source activate mirna

###########################################################################################
### Specify directories
###########################################################################################

IODIR=$1;shift
SCRIPTDIR='/diskmnt/Projects/Users/ssethura/cptac_mirna/Scripts/'
ANNDIR='/diskmnt/Projects/Users/ssethura/cptac_mirna/Annotation/'
REFFILE='/Projects/Users/qgao/Tools/refdata-GRCh38-2.1.0/fasta/genome.fa'

###########################################################################################
### Move to IO directory and start to loop over all ".fastq.gz files". Typically only one.
###########################################################################################

cd $IODIR

for t in *.fastq.gz
do
	f=${t%%.fastq.gz}""

###########################################################################################
### Adapter trimming
###########################################################################################

trimmomatic SE -threads 4 -trimlog ${f}_trimmomatic.log ${f}.fastq.gz ${f}.trimmed.fastq ILLUMINACLIP:${ANNDIR}/Adapters.fa:0:30:1 SLIDINGWINDOW:10:20 MINLEN:15 HEADCROP:4 TAILCROP:4 AVGQUAL:30

###########################################################################################
### QC report before and after trimming
###########################################################################################

fastqc ${f}.fastq.gz &
fastqc ${f}.trimmed.fastq &

###########################################################################################
### Alignment with BWA followed by splitting multiple alignments to individual lines in the SAM file
###########################################################################################

bwa aln -t 20 ${REFFILE} ${f}.trimmed.fastq > ${f}.sai 2> ${f}.bwa.log
bwa samse -n 10 ${REFFILE} ${f}.sai ${f}.trimmed.fastq > ${f}.sam
perl ${SCRIPTDIR}/xa2multi.pl ${f}.sam > temp.sam
mv temp.sam ${f}.sam

###########################################################################################
### Convert to BAM and BED and sort
###########################################################################################

samtools view -Sb ${f}.sam > ${f}.bam
bedtools bamtobed -i ${f}.bam > ${f}.bed
sort -k1,1 -k2,2n ${f}.bed > ${f}.sorted.bed

###########################################################################################
### Annotation using bedtools
###########################################################################################

### Annotate miRNAs

bedtools intersect -sorted -wao -f 8E-9 -a ${f}.sorted.bed -b ${ANNDIR}/mirbase.bed > ${f}_mirbase.bed
grep -v "ID" ${f}_mirbase.bed | cut -f1-6 > ${f}_unannot.bed
grep "ID" ${f}_mirbase.bed > temp 
mv temp ${f}_mirbase.bed

### Annotate tRNAs

bedtools intersect -sorted -wao -f 8E-9 -a ${f}_unannot.bed -b ${ANNDIR}/gencode_tRNAs.bed > ${f}_tRNAs.bed
grep -v "ID" ${f}_tRNAs.bed | cut -f1-6 > ${f}_unannot.bed
grep "ID" ${f}_tRNAs.bed > temp
mv temp ${f}_tRNAs.bed

### Annotate exons (includes other small RNAs such as rRNAs, snRNAs, snoRNAs, etc)

bedtools intersect -sorted -wao -f 2E-9 -a ${f}_unannot.bed -b ${ANNDIR}/gencode_std_exons.bed > ${f}_exons.bed
grep -v "ID" ${f}_exons.bed | cut -f1-6 > ${f}_unannot.bed
grep "ID" ${f}_exons.bed > temp
mv temp ${f}_exons.bed

###########################################################################################
### miRNA annotation prioritization
###########################################################################################

python ${SCRIPTDIR}/mirna_annotation.py --sample $f

grep "@" ${f}.sam > header.txt
cat ${f}_Annotated.sam >> header.txt
mv header.txt ${f}_Annotated.sam
sed "s/\t\n/\n/g" ${f}_Annotated.sam > temp.sam
sed "s/\t$//g" temp.sam > ${f}_Annotated.sam

###########################################################################################
### miRNA counting 
###########################################################################################

python ${SCRIPTDIR}/mirna_counting.py --sample $f --mirbase ${ANNDIR}/hsa.bed

###########################################################################################
### Close the loop 
###########################################################################################

done
