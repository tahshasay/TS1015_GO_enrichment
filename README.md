# TS1015_GO_enrichment

These data and scripts are supplement to:

Say & Degnan (in 1st review) Molecular and behavioural evidence that interdependent photo- and chemo-sensory systems regulate larval settlement in a marine sponge. Molecular Ecology. 

A preprint is available at: Say, T. E., & Degnan, S. M. (2019). Interdependent photo- and chemosensory systems regulate larval settlement in a marine sponge. bioRxiv, 519512. doi:10.1101/519512



#DEG were characterised by Gene Ontology (GO) enrichment using a Fishers exact test (gene is present or absent in the differentially expressed gene list) implemented using the GO_MWU package (Wright, Aglyamova, Meyer, & Matz, 2015); the package and instructions are available at https://github.com/z0on/GO_MWU and see additional supplementary data listed below.


#This folder contains the files and scripts used to perform the GO enrichment on the genes identified by the sPLS-DA.

INPUT_files
"Aq_all_V1.txt" - Aqu2.1 annotation file.

"splsda_PC1_ch3_SupFile3.7_1s.txt" - list of genes identified by the sPLS-DA on the 1st component. 

"gene2go.tab" - two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.

"GO_1&0s_splsda_PC1_ch3_SupFile3.7.txt" - correct format input file for analysis (includes DEG identified by the sPLS-DA (component 1). 

OUTPUT_files
GO enrichment plot is in Figure S4.

Reference
Wright, R. M., Aglyamova, G. V., Meyer, E., & Matz, M. V. (2015). Gene expression associated with white syndromes in a reef building coral, Acropora hyacinthus. BMC Genomics, 16(1), 371. doi:10.1186/s12864-015-1540-2
