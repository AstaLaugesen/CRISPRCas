#Making heat map of bacteria vs subtypes in the children:

#Loading phyloseq object of CRISPRs:
load("~/data/CRISPRs.RData")
#Loading bacterial data:
load("~/data/16S/F1y.RData")
#OTU table has 6846 taxa and 625 samples:
F1y

#Checking out the bacterial data further:
#head(sample_data(F1y))
#sample_variables(F1y)
#head(tax_table(F1y))

#agglomerate bacterial data by genus level
agglomerated_F1y <- tax_glom(F1y, taxrank=rank_names(F1y)[6])
#OTU table now has 669 taxa and still 625 samples
agglomerated_F1y

#most abundant bacteria in the sample (absolute abundance, not relative according to sequencing depth)
sort(taxa_sums(agglomerated_F1y), decreasing=TRUE)[0:10]
#Cannot conclude on this but it gives us an idea of what's probably most abundant



#For doing relative abundance, relatively to each row's sum:
install_github("JStokholm/Abundance")
library(abundance)

#taking relative abundance on bacteria table and CRISPR-Cas systems subtype table
genus_ab <- abundance(phylo_ob=F1y, level="Genus", id="Patient",relative_abun=TRUE, remove_collapsed_taxa=FALSE, select_taxa=NULL,select_level=NULL)
CRISPR_ab <- abundance(CRISPRs,level="Subtype",id = "ABCno",sample_id="SampleName")


#Sorting and cleaning up data:
#Ensuring we only have data in which ABC no. is in both tables
CRISPR_ab <- CRISPR_ab[CRISPR_ab$ABCno %in% genus_ab$Patient, ]
CRISPR_ab <- CRISPR_ab[!duplicated(CRISPR_ab$ABCno),]
CRISPR_ab <- CRISPR_ab[order(CRISPR_ab$ABCno),]

#doing the same for the bacterial data table
genus_ab <- genus_ab[genus_ab$Patient %in% CRISPR_ab$ABCno, ]
genus_ab <- genus_ab[order(genus_ab$Patient),]


#Taking out the data with the bacterias' abundance (first 3 columns being additional info)
abund_bac <- genus_ab[,3:ncol(genus_ab)]
#Ordering by median to find the most abundant bacteria
abund_bac <- abund_bac[,order(-apply(abund_bac, 2, median))]
#Selecting the 50 most abundant bacteria, and doing a log transformation
abund_bac <- abund_bac[,1:50] %>%  apply(2, function(x) log(x + min(x[x > 0])/2))

#removing _ from bacteria names for better layout in heatmap
colnames(abund_bac) <- gsub('^_','',colnames(abund_bac))
colnames(abund_bac) <- gsub('_',' ',colnames(abund_bac))


#Doing the same for the CRISPR data
abund_cris <- CRISPR_ab[,3:ncol(CRISPR_ab)]
#Ordering by median
abund_cris <- abund_cris[,order(-apply(abund_cris, 2, median))] 
#Select all system subtypes and do log transformation of the abundance
abund_cris <- abund_cris[,1:ncol(abund_cris)] %>%  apply(2, function(x) log(x + min(x[x > 0])/2))



library(pheatmap)

#Using Spearman's correlation to find the correlation between abundances of bacteria and CRISPR-Cas system subtypes
#And then making a heatmap
plot <- cor(abund_cris, abund_bac, method="spearman")%>% pheatmap(breaks = seq(-1,1,0.02), color = colorRampPalette(c("#000030", "darkblue", "lightblue", "white", "red", "darkred", "#300000"))(100), display_numbers = FALSE, fontsize_col = 14, fontsize_row = 14, legend = FALSE, annotation_legend = FALSE)
