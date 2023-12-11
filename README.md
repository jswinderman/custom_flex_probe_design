# custom_flex_probe_design
## Purpose
1. Import and subset references for custom probe targets.
2. Nominate probes and filter according to 10x custom probe guidelines.
3. Output custom probe order.

#### Author
Written by Jason Swinderman\
Reviewed by Aidan Winters

#### Date
Written November, 8th 2023\
Last revised: December, 4th 2023

#### 10x Custom Guide Design Rules
When designing custom probes for either singleplex or multiplex experiments, consider the following: 

• GC content should be between 44 − 72% for each 25 bp probe half. 

• Avoid homopolymer repeats. 

• Avoid overlap with annotated repeat or low complexity sequences. 

• If possible, design probes for coding regions of mRNA as opposed to untranslated regions. 

• The 25th nucleotide of the probe (3' most nucleotide of the LHS probe) must be a T. The opposing nucleotide in the target RNA must be an A. 

• Avoid common single nucleotide polymorphisms (SNPs) and potential mismatches at the ligation junction. Refer to the UCSC Genome Browser and the Single Nucleotide Polymorphism Database (dbSNP). If avoiding SNPs is not possible, SNPs and mismatches should be at least four bp away from the ligation junction.

• If probes can bind to sequences other than the target mRNA sequence, an off-target signal may be observed. To check for off-target homology, align the probe sequence to the reference transcriptome using the Basic Local Alignment Search Tool (BLAST). Matches to off-target genes should have at least five mismatches in at least one of the LHS or RHS probes to prevent efficient hybridization.