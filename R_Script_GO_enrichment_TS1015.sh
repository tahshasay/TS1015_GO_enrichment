
####################################################################################
# The GO enrichment analysis + script was published by Mikhail Matz
# https://github.com/z0on/GO_MWU
#
####################################################################################


#-----------------------------------------------------------------------------------------
# Instructions
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to iBPose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fraction of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu
#-----------------------------------------------------------------------------------------




#########################################################################################
# 1) Preparing files for input into the GO enrichment analysis
# 
# this script uses a list of genes of interest (module / DEG list) and annotates these with a 1 (important) or kme values (when using modules from WGCNA). Then compares this list of genes to a background list (genes in gemone) and assigns these genes with a 0 (therefore genes of interest are annotated with a 1 and all others are annotated with a 0). 

setwd("/opsin/u/tahsha/TSay_R_wd/Matz_GO_a00.01_AQlist/R_wd")
getwd()

#created above
allannot = read.delim("Aq_all_V1.txt", header = T, sep = "\t")

dim(allannot)
#18390 - all genes input into WGCNA

head(allannot)

#created r00.03
kme = read.delim("/opsin/u/tahsha/TSay_R_wd/Matz_GO_a00.01_AQlist/R_wd/splsda_PC1_ch3_SupFile3.7_1s.txt", header = F, sep = "\t")

dim(kme)
# 2414
head(kme)
#                V1          V2
# 1 Aqu2.1.00021_001 0.383047666
# 2 Aqu2.1.00153_001 0.296222804
# 3 Aqu2.1.00324_001 0.324178302




# merge 2 data tables (max)

# merge by the "Aq_model" field
# include all lines with all.x=T,all.y=T
splsda =merge(kme, allannot, by="V1", all.x=T,all.y=T)

write.table(as.data.frame(splsda),file="GO_1&0s_splsda_PC1_ch3_SupFile3.7.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')

# manually rm NA and  specify 0 for genes not included in the module and the kME value (number between 0 and 1, module membership score) for genes included in the module.

#in R wd
write.csv(as.data.frame(splsda),file="GO_1&0s_splsda_PC1_ch3_SupFile3.7.csv")

quit()
n

cd /opsin/u/tahsha/TSay_R_wd/

mkdir wd_Matz_GO_splsda_BP_20181218
mkdir wd_Matz_GO_splsda_MF_20181218
mkdir wd_Matz_GO_splsda_CC_20181218
#GO: B P (biological process), GO: M F (molecular function / GO: C C cellular component


#################################
# manually
# and remove first col
# open in textwrangler - replace NA with 0 (15976 replacements turq) save as _rmNAs
# save as .unix as filename_rmNSAs
#
#################################

# eg: goal
# Aqu2.1.00769_001	0.234416469
# Aqu2.1.00785_001	0.415432361
# Aqu2.1.00847_001	-0.122480935
# Aqu2.1.00867_001	0.321508861
# Aqu2.1.00870_001	0.282685696
# Aqu2.1.00880_001	0.17467574
# Aqu2.1.00890_001	0.235066643



##########################################################################################
#2) GO enrichment analysis
#
# This script was published by Mikhail Matz
# https://github.com/z0on/GO_MWU
#
# for this script you need to
#install.packages("ape")
# ‘/mnt/tBPdir/teBP/RtBPAzK1kW/downloaded_packages’

setwd("/opsin/u/tahsha/TSay_R_wd/wd_Matz_GO_splsda1_BP_20181219")

#laptop
#setwd("/Users/tahshasay/Documents/Scripts/R/wd_20181218/")
getwd()
input="GO_splsda_PC1_ch3_SupFile3.7_rmNAs.csv" # master_rmNAs.csv was renamed to splsda_rmNAs.csv on 20181218
#ExaBPle:
#input="blue_module.csv"
#
# my data
# TCONS_00005054	0.000219566
# TCONS_00005070	0.358856201
# TCONS_00005123	0.039197704
# Aqu2.1.00014_001	0
# Aqu2.1.00015_001	0
# Aqu2.1.00017_001	0
# Aqu2.1.00018_001	0

#input="heats.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
# GOI is largest file (probably all genes in the genome)
# wc -l ./*tab
# 44364 ./heats.csv


goAnnotations="gene2go.tab"
#goAnnotations="amil_defog_iso2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
# genes to go file = smaller than the GOI list
# 13381 ./amil_defog_iso2go.tab 
#
# manually removed: GO:P (biological process), GO:F (molecular .. in text wrangler
#read.table("matz_go.obo", sep = "\t") - error
goDatabase="go.obo" # download fromhttp://www.geneontology.org/GO.downloads.ontology.shtml


##################################################
#to change 
goDivision="BP" # either $BP$, or $BP$, or $BP$
##################################################


source("gomwu.functions.R")

library("ape")	


##################################################################
gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="/usr/bin/perl", 
	largest=0.1,  
	smallest=1,   # change to 1 because 1 TF can play iBP role
	clusterCutHeight=0.25)
	#,Module=TRUE,Alternative="g"
##################################################################


results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=0.001,#must be 0.001 for GO or 0.1 for other
	level1=1,   # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # originally 0.1
	level2=0.01,   # FDR cutoff to print in regular (not italic) font.
	level3=0.001,  # FDR cutoff to print in large bold font.
	txtsize=1.2,   # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5,# height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# manually rescale the plot so the tree matches the text 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

save(gomwuStats, goDatabase, input, goAnnotations, goDivision, results, file = "gomwuStats_BP_splsda_20181218.RData")

##########################################################################################
# 3) visualising the results


# to change
setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS0616_WGCNA_GO/b_enrichment_ouput/wd_Matz_GO_splsda2_BP_20181219")

# to change module v splsd
load(file = "gomwuStats_BP_splsda_20181218.RData")
load(file = "gomwuStats_CC_splsda_20181218.RData")
# for module:
# load(file = "gomwuStats_CC_turquoise_20181217.RData")

source("gomwu.functions.R")




library("ape")


#########################
# plot all GO categories
#########################

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=0.001,#must be 0.001 for GO or 0.1 for other
	level1=1,   # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # originally 0.1
	level2=0.01,   # FDR cutoff to print in regular (not italic) font.
	level3=0.001,  # FDR cutoff to print in large bold font.
	txtsize=0.5,   # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5,# height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)



##################
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=0.001,#must be 0.001 for GO or 0.1 for other - doesnt change anything
	level1=0.1,   # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # originally 0.1
	level2=0.01,   # FDR cutoff to print in regular (not italic) font.
	level3=0.001,  # FDR cutoff to print in large bold font.
	txtsize=1.2,   # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5,# height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)

quartz.save("splsda2_BP_20181219.pdf", type="pdf");


########################################################################
# text representation of results, with actual adjusted p-values
results
# color = colour on heatmap (not module colour)

write.table(results,"B23076_gomwuStats_BP_splsda2_20181222_results.txt",sep="\t", row.names=T, col.names=T) 
########################################################################


##################
# default settings
#
##################
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=0.001,#must be 0.001 for GO or 0.1 for other
	level1=0.1,   # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # originally 0.1
	level2=0.05,   # FDR cutoff to print in regular (not italic) font.
	level3=0.01,  # FDR cutoff to print in large bold font.
	txtsize=1.2,   # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5,# height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
quartz.save("splsda2_BP_20181219_matz_settings.pdf", type="pdf");

# manually rescale the plot so the tree matches the text 

########################################################################
# text representation of results, with actual adjusted p-values
results
# color = colour on heatmap (not module colour)

write.table(results,"B23076_gomwuStats_BP_splsda2_20181222_results_matz.txt",sep="\t", row.names=T, col.names=T) 
########################################################################


quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=0.001,#must be 0.001 for GO or 0.1 for other
	level1=0.05,   # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # originally 0.1
	level2=0.01,   # FDR cutoff to print in regular (not italic) font.
	level3=0.005,  # FDR cutoff to print in large bold font.
	txtsize=1.2,   # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5,# height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)


quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=0.001,#must be 0.001 for GO or 0.1 for other
	level1=0.5, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue. # originally 0.1
	level2=0.01, # FDR cutoff to print in regular (not italic) font.
	level3=0.001, # FDR cutoff to print in large bold font.
	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5, # height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
results
# color = colour on heatmap (not module colour)

write.table(results,"B23076_gomwuStats_BP_splsda2_20181222_results.txt",sep="\t", row.names=T, col.names=T) 


