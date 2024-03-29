CPTAC miRNA-Seq analysis
========================

### github: <https://github.com/ding-lab/CPTAC_miRNA>

#### Developed by Sunantha Sethuraman
#### Modified and maintained by Lijun Yao lijunyao@wustl.edu

------------------------------------------------------------------------

Processing description
----------------------

The raw data was made available as _realn.bam files.

------------------------------------------------------------------------

### Bam2bed:

Aligned bam was converted to sam files, which was then converted to a BED
file and sorted.

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
