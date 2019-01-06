CPTAC miRNA-Seq analysis
========================

### github: <https://github.com/ding-lab/CPTAC_miRNA>

#### Sunantha Sethuraman

------------------------------------------------------------------------

Processing description
----------------------

The raw data was made available as .fastq.gz files.

------------------------------------------------------------------------

### Annotation pre-processing:

Annotation information to be used in the pipeline were downloaded from
miRBase v22 and GENCODE v29. The downloaded GTFs were converted to BED
format files. For GENCODE, only the transcript variant labeled with the
'Basic' tag was used since it is the predominant transcript variant.
Annotations were limited to standard chromosomes Chr 1-22, X,Y and MT.
See annotation_prep.sh for the details. The annotations used in the
processing of CPTAC data (v1.0) is provided as Annotations.zip.

------------------------------------------------------------------------

### Adapter trimming:

The fastq.gz files were first trimmed using TRIMMOMATIC (Bolger, A. M.,
Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for
Illumina Sequence Data. Bioinformatics, btu170) with the adapter
sequences provided in Adapters.fa. The following constraints were used
during trimming: 1. Random 4 nucleotides at the start and end of the
read after adapter trimming were cropped; 2. Any read with average read
quality &lt; 30 was dropped; 3. Any read with average base quality &lt;
20 in any sliding window of 10 bases was dropped; 4. Reads shorter than
15 bases after processing were dropped.

------------------------------------------------------------------------

### Quality check:

FASTQC (<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>)
was used to quality check the reads before and after trimming to ensure
adapter removal.

------------------------------------------------------------------------

### Alignment:

The reads were aligned to the human genome (Reference GRCh38) using bwa
aln, allowing no mismatches. For any given read up to 10 alignments were
permitted. The aligned SAM file contains one line per read with the
information about multiple alignments reported in that line. For
downstream analysis purposes, it was required to have one alignment per
line and this was accomplished using the xa2multi.pl script. The SAM
file was then converted to a BED file and sorted.

------------------------------------------------------------------------

### Annotation:

The BED file was annotated using BEDTOOLS 'intersect' and several
preprocessed annotation files (see annotation pre-processing above). A
custom python script (mirna_annotation.py) was used to handle multiple
annotations and multi-mapping of the reads to obtain one unique mapping
per read and one unique annotation per mapping.

------------------------------------------------------------------------

### Multiple annotations were resolved in the order of priority:

'miRNA', 'miRNA_primary_transcript', 'tRNA', 'rRNA', 'sRNA', 'snRNA',
'snoRNA', 'scaRNA', 'scRNA', 'ribozyme', 'vaultRNA', 'Mt_tRNA',
'Mt_rRNA', 'protein_coding','3prime_overlapping_ncRNA',
'bidirectional_promoter_lncRNA', 'lincRNA','macro_lncRNA',
'misc_RNA', 'non_coding','antisense', 'sense_intronic',
'sense_overlapping', 'IG_C_gene', 'IG_D_gene', 'IG_J_gene',
'IG_V_gene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene',
'TR_V_gene', 'IG_C_pseudogene','IG_J_pseudogene',
'IG_V_pseudogene', 'IG_pseudogene', 'TR_J_pseudogene',
'TR_V_pseudogene', 'rRNA_pseudogene', 'processed_transcript',
'pseudogene', 'polymorphic_pseudogene', 'processed_pseudogene',
'transcribed_processed_pseudogene',
'transcribed_unitary_pseudogene',
'transcribed_unprocessed_pseudogene',
'translated_processed_pseudogene', 'unitary_pseudogene',
'unprocessed_pseudogene', 'TEC','unannotated'

Since mature miRNAs are processed from precursor miRNAs (referred above
as 'miRNA_primary_transcript'), most of the reads mapping to mature
miRNAs also map to precursor miRNAs and vice-versa. The priorities were
assigned by overlap length, allowing a 3 nt additional advantage to
mature miRNAs. This advantage allow to correct for processing errors
(biological and sequencing based) and the number 3 was empirically
determined.

Multimapping was resolved based on the same order of priority as above.
Additionally, multi-maps within the same annotation group were assigned
randomly to one location.

------------------------------------------------------------------------

### Couting:

For each mature miRNA and precursor miRNA, the number of reads were
counted. The raw counts were then converted to TPM (Transcripts per
million) values using (raw_counts)*10e6/(total_counts), where
total_counts = total number of reads aligned to all miRNAs (mature +
precursor). It is important to note that TPM was not calculated as a
function of library size. The counts were summed for each precursor
miRNA (precursor + isoforms) and is reported as one file.

------------------------------------------------------------------------

### Note on miRNA orthologs:

Some mature miRNAs are orthologs of each other, meaning that the mature
miRNA sequence is the same, even if they arise from two distinct
precursors. For example, MIR26A1 and MIR26A2 are two genes, and the
miR-26a-5p that are derived from both these genes are identical in their
sequence and hence presumably their functions. miRBase distinguishes
these mature orthologs using their IDs. For example:

Percursor miR-26a-1 uses the ID MI0000083 Precursor miR-26a-2 uses the
ID MI0000750 Mature miR-26a-5p (derived from miR-26a-1) uses the ID
MIMAT0000082_1 Mature miR-26a-5p (derived from miR-26a-2) uses the ID
MIMAT0000082

Thus miRNA IDs serve as a better means of tracking orthologs. In this
analysis, all tables are designed to have one unique row per miRNA ID.
If users are interested in miR-26a-5p irrespective of its gene of
origin, it would be necessary to sum the values from the rows
corresponding to MIMAT0000082 and MIMAT0000082_1.

It is also important to note that during sequencing, these two mature
miRNAs (MIMAT0000082 and MIMAT0000082_1) would be practically
indistingushable. The pipeline is desgined to randomly assign such reads
to one location. So one would expect to have roughly similar counts for
miRNA orthologs. While MIR26A is one example, there exist other
orthologs and handling those are left to discretion of the user.

------------------------------------------------------------------------

### Output:

For each sample, three output files are provided with the following
suffixes: ".mature.txt", ".precursor.txt" and ".total.txt" which contain
information on mature miRNA counts, precursor miRNA counts and total
miRNA counts, respectively.

------------------------------------------------------------------------

### Sample nomenclature

The samples follow the naming system "[SubjectID].[T or
A].type.txt", where T or A specifies whether it is a tumor or a
tumor_normal sample (for example: C3N-01521.T.mature.txt). Type refers
to mature, precursor or total.
