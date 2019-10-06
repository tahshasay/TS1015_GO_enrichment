# TS1015_GO_enrichment

#These data and scripts are supplement to Say, T. E., & Degnan, S. M. (2019). Interdependent photo- and chemosensory systems regulate larval settlement in a marine sponge. bioRxiv, 519512. doi:10.1101/519512



#This folder contains the files and scripts used to perform the GO enrichment on the genes identified by the sPLS-DA.

INPUT_files
"Aq_all_V1.txt" - Aqu2.1 annotation file.

"splsda_PC1_ch3_SupFile3.7_1s.txt" - list of genes identified by the sPLS-DA on the 1st component. 

"gene2go.tab" - two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.

"GO_1&0s_splsda_PC1_ch3_SupFile3.7.txt" - correct format input file for analysis (includes DEG identified by the sPLS-DA (component 1). 

OUTPUT_files
GO enrichment plot is in Figure S4.
