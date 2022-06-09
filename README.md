## This project has 2 goals:

1. Speed up (improve) the old Perl scripts;

2. Translate all the Perl scripts into Python (3.6 or higher).

## Important bug fixes:

https://github.com/f-i/pySNPAAMapper/commit/a7dd5803d9165380fa3cf68c3997083056ed12f6 (Commit on Aug 19, 2020)

https://github.com/f-i/pySNPAAMapper/commit/1fda3b1b856efddfedf1c5fa91e0d61f24904c11 (Commit on May 28, 2022)

https://github.com/f-i/pySNPAAMapper/commit/5f76ca4ae3de9b428e78c609557d6174b07496f8 (Commit on Jun 1, 2022)

---
- Authors: Bai, Y & J. Cavalcoli
- SNPAAMapper Version 2.0 (Copyright 2013)
- Python Translator: Chenjian Fu
---

## Program descriptions:
SNPAAMapper is a downstream variant annotation program that can effectively classify variants by region (e.g. exon, intron, etc), predict amino acid change type (e.g. synonymous, non-synonymous mutation, etc), and prioritize mutation effects (e.g. CDS versus 5’UTR, etc). pySNPAAMapper is a Python version of SNPAAMapper.

## Major features:
1. The pipeline accepts the VCF input file in tab-delimited format and processes the vcf input file containing all cases (G5, lowFreq, and novel)
2. The variant mapping step has the option of letting users select whether they want to report the bp distance between each identified intron variant and its nearby exon
3. The pipeline can deal with VCF files called by different SAMTools versions (0.1.18 and older ones) and also offers flexibility in dealing with vcf input files generated using SAMTools with two or three samples
4. The spreadsheet result file contains full protein sequences for both ref and alt alleles, which makes it easier for downstream protein structure/function analysis tools to take


## Instructions:

Please dump all files in the same directory on Unix or Mac machines. The user can simply type
```
./run_SNPAAMapper.sh config.txt
```
as a commandline, or run the following steps in a sequential order (Note: the first two steps were compiled for human hg19 genome and output files were generated already):

### 1) Generate exon annotation file
```
perl Algorithm_generating_annotation_exon.pl ChrAll_knownGene.txt
```
or
```
python algorithm_generating_annotation_exon.py ChrAll_knownGene.txt
```

### 2) Process exon annotation files and generate feature start and gene mapping files
```
perl Algorithm_preprocessing_exon_annotation_RR.pl ChrAll_knownGene.txt.exons
```
or
```
python algorithm_preprocessing_exon_annotation_rr.py ChrAll_knownGene.txt.exons
```

### 3) Classify variants by regions (CDS, Upstream, Downstream, Intron, UTRs...)
This Perl script is extremely slow. The Python script is much faster.
```
perl Algorithm_mapping_variants_reporting_class_intronLocation_updown.pl ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format
```
or
```
python algorithm_mapping_variants_reporting_class_intronlocation_updown.py ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format
```
(e.g. use 007_crop.vcf as the VCF_input_file_in_tab_delimited_format)
If IntronExon_boundary_in_bp is considered:
```
perl Algorithm_mapping_variants_reporting_class_intronLocation_updown.pl ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format IntronExon_boundary_in_bp
```
or
```
python algorithm_mapping_variants_reporting_class_intronlocation_updown.py ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format IntronExon_boundary_in_bp
```
(e.g. use 007_crop.vcf and 6 as the VCF_input_file_in_tab_delimited_format and IntronExon_boundary_in_bp)

### 4) Predict amino acid change type
```
perl Algorithm_predicting_full_AA_change_samtools_updown.pl VCF_input_file_in_tab_delimited_format.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.out.txt
```
or
```
python algorithm_predicting_full_aa_change_samtools_updown.py VCF_input_file_in_tab_delimited_format.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.out.txt
```
(e.g. use 007_crop.vcf.append and 007_crop.vcf.out.txt as the VCF_input_file_in_tab_delimited_format.append and VCF_input_file_in_tab_delimited_format.out.txt)

Note: The hg19_CDSIntronWithSign.txt.out file's size is huge and usually around
3 GB, which exceeds GitHub's file size limit of 100.00 MB. So it is not in this
repository. Please download it in "Source Code" from:
http://isu.indstate.edu/ybai2/SNPAAMapper2/index.html

### 5) Prioritize mutation effects
```
perl Algorithm_prioritizing_mutation_headerTop_updown.pl VCF_input_file_in_tab_delimited_format.append.out.txt
```
or
```
python algorithm_prioritizing_mutation_headertop_updown.py VCF_input_file_in_tab_delimited_format.append.out.txt
```
(e.g. use 007_crop.vcf.append.out.txt as the VCF_input_file_in_tab_delimited_format.append.out.txt)


> The final output file is \*.append.out.txt.prioritized_out.



## Terminology

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

CDS (coding sequence), the coding region of a gene, is the portion of a gene's
DNA or RNA that codes for protein.
