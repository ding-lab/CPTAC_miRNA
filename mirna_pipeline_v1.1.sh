###########################################################################################
### Activate conda environment (ignore if all the necessary softwares are installed in the local environment)
###########################################################################################

source activate mirna

###########################################################################################
### Obtain input and specify directories
### This script take Unaln bams as inputs #Lijun modified from Sunantha's pipeline
###########################################################################################

INFILE=$1;shift
OUTDIR=$1;shift
SCRIPTDIR='/diskmnt/Projects/Users/ssethura/cptac_mirna/Scripts/'
ANNDIR='/diskmnt/Projects/Users/ssethura/cptac_mirna/Annotation/'
#REFFILE='/Projects/Users/qgao/Tools/refdata-GRCh38-2.1.0/fasta/genome.fa' ## Denali
REFFILE='/diskmnt/Datasets/Reference/CellRanger/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa' ## Katmai

###########################################################################################
### Create output and temp directories. Clean up file name to get sample name.
###########################################################################################

INDIR=$(dirname "$INFILE")
fname=$(basename "$INFILE")
dname=$(basename $(dirname "$INFILE"))

WRITEDIR=${OUTDIR}/${dname}_out
TMPDIR=${OUTDIR}/${dname}_tmp

mkdir ${WRITEDIR}
mkdir ${TMPDIR}

f=${fname%%.unaln.bam}""

###########################################################################################
### Move to TMP directory for all operations
###########################################################################################

cd $TMPDIR

###########################################################################################
### Ualn BAM to FASTQ
###########################################################################################
samtools bam2fq ${INDIR}/${f}.unaln.bam | gzip > ${f}.fastq.gz 2> ${f}.bam2fastq.log

###########################################################################################
### Adapter trimming
###########################################################################################

trimmomatic SE -threads 4 -trimlog ${f}_trimmomatic.log ${f}.fastq.gz ${f}.trimmed.fastq ILLUMINACLIP:${ANNDIR}/Adapters.fa:0:30:1 SLIDINGWINDOW:10:20 MINLEN:15 HEADCROP:4 TAILCROP:4 AVGQUAL:30

###########################################################################################
### QC report before and after trimming
###########################################################################################

fastqc ${f}.fastq.gz -o ${TMPDIR}
fastqc ${f}.trimmed.fastq -o ${TMPDIR}

###########################################################################################
### Alignment with BWA
###########################################################################################

bwa aln -t 20 ${REFFILE} ${f}.trimmed.fastq > ${f}.sai 2> ${f}_bwa.log
bwa samse -n 10 ${REFFILE} ${f}.sai ${f}.trimmed.fastq > ${f}.sam

###########################################################################################
### Write the number of reads in FASTQs and SAM to the report
###########################################################################################

echo "Number of reads in FASTQ" > ${f}_report.txt
zcat ${f}.fastq.gz | grep "@" | wc -l >> ${f}_report.txt
echo "Number of reads in trimmed FASTQ" >> ${f}_report.txt
grep "@" ${f}.trimmed.fastq | wc -l >> ${f}_report.txt
echo "Number of reads in SAM" >> ${f}_report.txt
grep -v "@" ${f}.sam | wc -l >> ${f}_report.txt
echo "Number of reads with multiple alignments in SAM" >> ${f}_report.txt
grep "XA:" ${f}.sam | wc -l >> ${f}_report.txt
echo "Number of reads with unique alignments in SAM" >> ${f}_report.txt
grep "NM:" ${f}.sam | grep -v "XA:" | wc -l >> ${f}_report.txt
echo "Number of unaligned reads in SAM" >> ${f}_report.txt
grep -v "NM:" ${f}.sam | grep -v "@" | wc -l >> ${f}_report.txt

###########################################################################################
### Split reads with multiple alignments into individual lines
###########################################################################################

perl ${SCRIPTDIR}/xa2multi.pl ${f}.sam > temp.sam
mv temp.sam ${f}.sam

###########################################################################################
### Convert to BAM and BED and sort
###########################################################################################

samtools view -Sb ${f}.sam > ${f}.bam
bedtools bamtobed -i ${f}.bam > ${f}.bed
sort -T ${TMPDIR} -k1,1 -k2,2n ${f}.bed > ${f}.sorted.bed

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

python ${SCRIPTDIR}/mirna_annotation.py --sample $f 2> ${f}_annotation.log

grep "@" ${f}.sam > header.txt
cat ${f}_Annotated.sam >> header.txt
mv header.txt ${f}_Annotated.sam
sed "s/\t\n/\n/g" ${f}_Annotated.sam > temp.sam
sed "s/\t$//g" temp.sam > ${f}_Annotated.sam

###########################################################################################
### miRNA counting 
###########################################################################################

python ${SCRIPTDIR}/mirna_counting.py --sample $f --mirbase ${ANNDIR}/mirbase.bed 2> ${f}_counting.log

###########################################################################################
### Count number of reads in annotated files and add to the report
###########################################################################################

#echo "Number of lines in annotated SAM (multiple alignments and multiple annotations each make a new line)"  >> ${f}_report.txt
#grep -v "@" ${f}_Annotated.sam | wc -l >> ${f}_report.txt
echo "Number of lines in annotations (one unique annotation per read)" >> ${f}_report.txt
grep -v "read" ${f}_Annotated.txt | wc -l >> ${f}_report.txt
echo "Number of mature miRNAs" >> ${f}_report.txt
grep -v "read" ${f}_mature.txt | wc -l >> ${f}_report.txt
echo "Number of precursor miRNAs" >> ${f}_report.txt
grep -v "read" ${f}_precursor.txt | wc -l >> ${f}_report.txt
echo "Number of total miRNAs" >> ${f}_report.txt
grep -v "read" ${f}_miRNA_expression.txt | wc -l >> ${f}_report.txt

###########################################################################################
### Final number check
###########################################################################################

#trimmed = $(grep "@" ${f}.trimmed.fastq | wc -l)
#annotated = $(grep -v "read" ${f}_Annotated.txt | wc -l)

#sam = $(grep -v "@" ${f}.sam | wc -l)
#asam = $(grep -v "@" ${f}_Annotated.sam | wc -l)

###########################################################################################
### Move output files to WRITEDIR and delete TMPDIR if read numbers tally
###########################################################################################

#if [ "$trimmed" == "$annotated" ]
#then
#if [ "$sam" == "$asam" ]
#then

mv -t ${WRITEDIR} ${f}_Annotated.txt ${f}_mature.txt ${f}_precursor.txt ${f}_miRNA_expression.txt {f}*fastqc.html ${f}_report.txt # ${f}_Annotated.sam

#rm -r ${TMPDIR}

#fi
#fi
