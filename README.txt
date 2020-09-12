########################################
Authors: Bai, Y & J. Cavalcoli
SNPAAMapper Version 2.0 (Copyright 2013)
Python Translator: Chenjian Fu (with 2 Goals here: 1 Speed up old Perl scripts; 2 Translate Perl into Python)
########################################

Program descriptions:
SNPAAMapper is a downstream variant annotation program that can effectively classify variants by region (e.g. exon, intron, etc), predict amino acid change type (e.g. synonymous, non-synonymous mutation, etc), and prioritize mutation effects (e.g. CDS versus 5’UTR, etc).

Major features:
a) The pipeline accepts the VCF input file in tab-delimited format and processes the vcf input file containing all cases (G5, lowFreq, and novel)
b) The variant mapping step has the option of letting users select whether they want to report the bp distance between each identified intron variant and its nearby exon
c) The pipeline can deal with VCF files called by different SAMTools versions (0.1.18 and older ones) and also offers flexibility in dealing with vcf input files generated using SAMTools with two or three samples
d) The spreadsheet result file contains full protein sequences for both ref and alt alleles, which makes it easier for downstream protein structure/function analysis tools to take


Instructions: Please dump all files in the same directory on Unix or Mac machines. The user can simply type
./run_SNPAAMapper.sh config.txt
or run the following steps in a sequential order (Note: the first two steps were compiled for human hg19 genome and output files were generated already):

# Generate exon annotation file
1) perl Algorithm_generating_annotation_exon.pl ChrAll_knownGene.txt
or
 python Algorithm_generating_annotation_exon.py ChrAll_knownGene.txt

# Process exon annotation files and generate feature start and gene mapping files
2) perl Algorithm_preprocessing_exon_annotation_RR.pl ChrAll_knownGene.txt.exons
or
 python Algorithm_preprocessing_exon_annotation_RR.py ChrAll_knownGene.txt.exons

# Classify variants by regions (CDS, Upstream, Downstream, Intron, UTRs...)
3) perl(/python) Algorithm_mapping_variants_reporting_class_intronLocation_updown.pl(/py) ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format.txt
(e.g. perl(/python) Algorithm_mapping_variants_reporting_class_intronLocation_updown.pl(/py) ChrAll_knownGene.txt.exons 007_crop.vcf)
OR
perl(/python) Algorithm_mapping_variants_reporting_class_intronLocation_updown.pl(/py) ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format.txt IntronExon_boundary_in_bp
(e.g. perl(/python) Algorithm_mapping_variants_reporting_class_intronLocation_updown.pl(/py) ChrAll_knownGene.txt.exons 007_crop.vcf 6)

# Predict amino acid change type
4) perl Algorithm_predicting_full_AA_change_samtools_updown.pl VCF_input_file_in_tab_delimited_format.txt.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.txt.out.txt
(perl Algorithm_predicting_full_AA_change_samtools_updown.pl 007_crop.vcf.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt >007_crop.vcf.out.txt)

# Prioritize mutation effects
5) perl Algorithm_prioritizing_mutation_headerTop_updown.pl VCF_input_file_in_tab_delimited_format.txt.append.out.txt
(perl Algorithm_prioritizing_mutation_headerTop_updown.pl 007_crop.vcf.append.out.txt)


The final output file is *.append.out.txt.prioritized_out.



Terminology:

In molecular genetics, an untranslated region (or UTR) refers to either of two
sections, one on each side of a coding sequence on a strand of mRNA. If it is
found on the 5' side, it is called the 5' UTR (or leader sequence), or if it is
found on the 3' side, it is called the 3' UTR (or trailer sequence).

By convention, upstream and downstream relate to the 5' to 3' direction
respectively in which RNA transcription takes place. Upstream is toward the 5'
end of the RNA molecule and downstream is toward the 3' end. When considering
double-stranded DNA, upstream is toward the 5' end of the coding strand for the
gene in question and downstream is toward the 3' end.

A chromosome (Ran Se Ti) is a DNA molecule with part or all of the genetic
material of an organism.

UTR3Start UTR5Start
UTR3End UTR5End
downstreamStart upstreamStart
downstreamEnd upstreamEnd
down_flank up_flank

Single nucleotide polymorphisms, frequently called SNPs (pronounced “snips”),
are the most common type of genetic variation among people.
