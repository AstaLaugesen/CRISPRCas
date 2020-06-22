#loading libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(readxl)
library(ggplot2)


#reading in files
setwd("~/data")
otu_mat<- t(read.table("crispr/S_allcrisprsystems_abundance_crisprcas_and_casorphans.txt", header=TRUE, sep="\t", row.names = 1))
tax_mat<- read_excel("crispr/taxonomy.xlsx", sheet = "Data")
fishoil <- read_excel("clinical/ABC0058_Fishoil.xlsx", sheet="Data")
sampleIDs <-read.table("clinical/samples.abc.tab", header=FALSE, sep="\t")
Dvitamin <-read_excel("clinical/ABC0059_D-vitamin.xlsx", sheet="Data")
delivery <-read_excel("clinical/ABC0016_Delivery.xlsx", sheet="Data")
animals <-read_excel("clinical/ABC0095_Furred_animals_1yr.xlsx", sheet="Data")
asthma <- read_excel("clinical/J45_cox_cross_180612.xlsx", sheet="Sheet 1")
guts <-read_excel("clinical/gut_scores_1y.xlsx", sheet="Sheet 1")
mother <-read.table("clinical/asthma.mother.tsv", header=TRUE, sep="\t")
antibiotics <-read.table("clinical/ab.tsv", header=TRUE, sep="\t")
ruralurban <-read.table("clinical/abc.ruralUrban.tsv", header=TRUE, sep="\t")

#checking rows and cols fit in the otu table
#row.names(otu_mat)
#colnames(otu_mat)

#defining row names in the taxonomy table
row.names(tax_mat) <- tax_mat$ID
tax_mat <- tax_mat %>% select (-ID) 

#removing the patient no. that wasn't sequenced (no sample number)
sampleIDs <- sampleIDs[-c(403),]
#adding the column with the S[...] identifiers and adding column titles
#(Done because samples can't start with a number)
sampleIDs$V0 <- paste("S", sampleIDs$V1, sep="")
colnames(sampleIDs) <- c("SampleName", "ABCno", "NewNames")

#Changing column names and removing comments/other less relevant variables
fishoil <- fishoil %>% select (-Comments)
colnames(fishoil) <- c("ABCno", "Supplement_oil")
Dvitamin <- Dvitamin %>% select (-Comments)
colnames(Dvitamin) <- c("ABCno", "Supplement_vitamin")
colnames(delivery) <- c("ABCno", "Birthdate", "Delivery","Sectiotype")
animals <- animals %>% select (-NUMBEROFCATS) %>% select (-NUMBEROFDOGS) %>% select (-FURREDANIMALOTHERSPEC) %>% select (-COMMENTS)
colnames(animals) <- c("ABCno", "Furred_animal_days", "Days_with_cat","Days_with_dog","Days_with_other_animals")
asthma <- select(asthma, abcno, j45_5yr_cross, j45_5yr_ever)
colnames(asthma) <- c("ABCno", "j45_5yr_cross", "j45_5yr_ever")
guts <- guts %>% select (-j45_cross_pls_1y_all)
colnames(guts) <- c("ABCno", "pamcluster", "maz")
colnames(mother) <- c("ABCno", "Asthmatic_mother")
colnames(antibiotics) <- c("ABCno", "Antibiotics_birth_child","Antibiotics_birth_mother", "Antibiotics_ever_1yr_child")
colnames(ruralurban) <- c("ABCno", "RuralUrbanStatus")

#left join the samples with ABC numbers
samples_df<-merge(x=sampleIDs,y=fishoil,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=Dvitamin,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=delivery,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=animals,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=asthma,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=guts,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=mother,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=antibiotics,by="ABCno",all.x=TRUE)
samples_df<-merge(x=samples_df,y=ruralurban,by="ABCno",all.x=TRUE)

#defining rows in sample tables
row.names(samples_df) <- samples_df$NewNames
samples_df <- samples_df %>% select (-NewNames)


#And now putting it all together:

#making OTU and tax table into matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform into phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAMPLES = sample_data(samples_df)

CRISPRs_unfiltered <- phyloseq(OTU, TAX, SAMPLES)

#Checking if there are samples in which no CRISPR-Cas systems were found with CCTyper (we already know there are)
sample_sums(CRISPRs_unfiltered)== 0
#We get it confirmed that sample 10 and 642 are missing CRISPR-Cas systems, so these are filtered out
CRISPRs <- subset_samples(CRISPRs_unfiltered, sample_names(CRISPRs_unfiltered) != "S18097D-02-10" & sample_names(CRISPRs_unfiltered) != "S18097D-02-642")

#Lastly saving to avoid having to run it through every time
save(CRISPRs, file = "CRISPRs.RData", compress = T)









library(devtools)
install_github("JStokholm/rabuplot")
library(rabuplot)

#Univariate: 


#Overview:
rabuplot(CRISPRs, "ABCno", "Type", main="Average relative abundance for all samples", bar_chart = TRUE, bar_chart_stacked = TRUE)
ggsave("plots/overview.pdf", width = 20, height =10 )

#Supplements:
##oil
rabuplot(CRISPRs, "Supplement_oil", "Subtype", main="CRISPR-Cas system subtypes in children with mother receiving fishoil or placebo", bar_chart = TRUE, bar_chart_stacked = TRUE)
#ggsave("plots/fishoil_barplot_subtypes.pdf", width = 11, height =5 )  
rabuplot(CRISPRs, "Supplement_oil", "Subtype", main="CRISPR-Cas system subtypes in children with mother receiving fishoil or placebo", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/fishoil_barplot_subtypes.pdf", width = 11, height =5 ) 
#I-B at 0.024

##vitamin
rabuplot(CRISPRs, "Supplement_vitamin", "Subtype", main="CRISPR-Cas system subtypes in children with mother receiving vitaminD or placebo", bar_chart = TRUE, bar_chart_stacked = TRUE)
rabuplot(CRISPRs, "Supplement_vitamin", "Subtype", main="CRISPR-Cas system subtypes in children with mother receiving vitaminD or placebo", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/vitamin_barplot_subtypes.pdf", width = 11, height =5 )
#III-D at 0.036

#Asthma
rabuplot(CRISPRs, "j45_5yr_ever", "Subtype", main="CRISPR-Cas system subtypes for children who have had asthma at age 5", bar_chart = TRUE, bar_chart_stacked = TRUE)
ggsave("plots/asthmaever_barplot_stacked_subtypes.pdf", width = 11, height =5 ) 
rabuplot(CRISPRs, "j45_5yr_cross", "Subtype", main="CRISPR-Cas system subtypes for children who had asthma at age 5", bar_chart = TRUE, bar_chart_stacked = TRUE)
ggsave("plots/asthmacross_barplot_stacked_subtypes.pdf", width = 11, height =5 ) 
rabuplot(CRISPRs, "j45_5yr_ever", "Subtype", main="CRISPR-Cas system subtypes for children diagnosed with asthma at any time between age 0-5", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/asthmaever_barplot_subtypes.pdf", width = 11, height =5 ) 
#VI-D at 0.035, II-B at 0.034
rabuplot(CRISPRs, "j45_5yr_cross", "Subtype", main="CRISPR-Cas system subtypes for children diagnosed with asthma at age 5", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/asthmacross_barplot_subtypes.pdf", width = 11, height =5 )
#none

#Delivery
rabuplot(CRISPRs, "Delivery", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
rabuplot(CRISPRs, "Delivery", "Subtype", main="CRISPR-Cas systems in children born with different methods", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/delivery_barplot_subtypes.pdf", width = 11, height =5 )
#II-A at 0.042

#Pamcluster - good seque into bacteria introduction
rabuplot(CRISPRs, "pamcluster", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
ggsave("plots/pamclustering_barplotStacked_subtypes.pdf", width = 11, height =5)
rabuplot(CRISPRs, "pamcluster", "Subtype", main="CRISPR-Cas system subtypes for children with different microbial community types (PAM clusters)", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/pamclustering_barplot_subtypes.pdf", width = 11, height =5)
#I-C w/ 0.042, II-C w/ <0.001, I-E w/ 0.039, III-A w/ 0.010, V-A w/ <0.001

#Asthmatic mother
rabuplot(CRISPRs, "Asthmatic_mother", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
rabuplot(CRISPRs, "Asthmatic_mother", "Subtype", main="CRISPR-Cas system subtypes for children with asthmatic mothers", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/asthmaticmother_barplot_subtypes.pdf", width = 11, height =5) 
#II-A with 0.011 should be pretty significant

#Antibiotics:
##child, at birth:
rabuplot(CRISPRs, "Antibiotics_birth_child", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
rabuplot(CRISPRs, "Antibiotics_birth_child", "Subtype", main="CRISPR-Cas system subtypes in children given antibiotics at birth", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/antibioticschildBirth_barplot_subtypes.pdf", width = 11, height =5)
#II-B at 0.030
##mother, at birth of child:
rabuplot(CRISPRs, "Antibiotics_birth_mother", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
rabuplot(CRISPRs, "Antibiotics_birth_mother", "Subtype", main="CRISPR-Cas system subtypes in children with mother given antibiotics during labor", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/antibioticsmother_barplot_subtypes.pdf", width = 11, height =5)
#VI-B at 0.019
##child, at any point up to 1yr:
rabuplot(CRISPRs, "Antibiotics_ever_1yr_child", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
rabuplot(CRISPRs, "Antibiotics_ever_1yr_child", "Subtype", main="CRISPR-Cas system subtypes in children given antibiotics between age 0-1", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/antibioticchildever1yr_barplot_subtypes.pdf", width = 11, height =5)
#II-A seems significant at 0.003, V-A seems significant at 0.014 

#Rural/Urban status
rabuplot(CRISPRs, "RuralUrbanStatus", "Subtype", main="", bar_chart = TRUE, bar_chart_stacked = TRUE)
ggsave("plots/RuralUrban_barplot_stacked_subtypes.pdf", width = 11, height =5) 
rabuplot(CRISPRs, "RuralUrbanStatus", "Subtype", main="CRISPR-Cas system subtypes for children with different home placement", bar_chart = TRUE, bar_chart_stacked = FALSE)
ggsave("plots/RuralUrban_barplot_subtypes.pdf", width = 11, height =5) 
#II-B with 0.035 



#Multivariant testing:


#prepping for adonis:
#General distance matrix:
CRISPRs.d = distance(CRISPRs, method = "bray")
#For oil:
index <- !is.na(get_variable(CRISPRs,"Supplement_oil"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Supplement_oil")[index])
#Pr(>F)=0.651

#For vitamin:
index <- !is.na(get_variable(CRISPRs,"Supplement_vitamin"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Supplement_vitamin")[index])
#Pr(>F)=0.772

#For delivery:
index <- !is.na(get_variable(CRISPRs,"Delivery"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Delivery")[index])
#Pr(>F)=0.027

#for asthma
#ever:
index <- !is.na(get_variable(CRISPRs,"j45_5yr_ever"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "j45_5yr_ever")[index])
#Pr(>F)=0.425
#cross:
index <- !is.na(get_variable(CRISPRs,"j45_5yr_cross"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "j45_5yr_cross")[index])
#Pr(>F)=0.099

#for pamcluster
index <- !is.na(get_variable(CRISPRs,"pamcluster"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "pamcluster")[index])
#Pr(>F)=0.001

#for asthmatic mother
index <- !is.na(get_variable(CRISPRs,"Asthmatic_mother"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Asthmatic_mother")[index])
#Pr(>F)=0.218


#for antibiotics

#child, at birth:
index <- !is.na(get_variable(CRISPRs,"Antibiotics_birth_child"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Antibiotics_birth_child")[index])
#Pr(>F)=0.021
#child, at birth (only normal birth):
index <- !is.na(get_variable(CRISPRs,"Antibiotics_birth_child")) & get_variable(CRISPRs, "Delivery") =="Normal"
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Antibiotics_birth_child")[index])
#Pr(>F)=0.004
#mother, at birth of child:
index <- !is.na(get_variable(CRISPRs,"Antibiotics_birth_mother"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Antibiotics_birth_mother")[index])
#Pr(>F)=0.001
#mother, at birth of child (only normal birth, since every c-section mother gets antibiotica):
index <- !is.na(get_variable(CRISPRs,"Antibiotics_birth_mother")) & get_variable(CRISPRs, "Delivery") =="Normal"
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Antibiotics_birth_mother")[index])
#Pr(>F)=0.003
#Child, at any time up until 1yr:
index <- !is.na(get_variable(CRISPRs,"Antibiotics_ever_1yr_child"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Antibiotics_ever_1yr_child")[index])
#Pr(>F)=0.012

#Furred animals
index <- !is.na(get_variable(CRISPRs,"Furred_animal_days"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "Furred_animal_days")[index])
#Pr(>F)=0.971

#Rural/Urban status
index <- !is.na(get_variable(CRISPRs,"RuralUrbanStatus"))
adonis2(CRISPRs.d %>% as.matrix %>% .[index,index] %>% as.dist ~ get_variable(CRISPRs, "RuralUrbanStatus")[index])
#Pr(>F)=0.774
#rural/urban status usually always has significance when it comes to the bacteria. Interesting that it doesn't here


#ordination of plot PCoA, bray-curtis
#Antibiotics administered to mother
antibiotics_mother_filter <- subset_samples(CRISPRs, !is.na(Antibiotics_birth_mother)) 
plot_ordination(antibiotics_mother_filter,ordination=ordinate(antibiotics_mother_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(Antibiotics_birth_mother))) +
  stat_ellipse(aes(color=factor(Antibiotics_birth_mother), fill=factor(Antibiotics_birth_mother)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="Antibiotics administered to mother at birth (P<0.001)     ",breaks=c("1","0"), labels=c("Yes", "No")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=16)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.6%)", x="PC1  (23.1%)")
ggsave("plots/antibioticsmother_PCoA_subtypes.pdf", width = 11, height =6)


#Antibiotics administered to mother, only natural births
antibiotics_mother_filter2 <- subset_samples(CRISPRs, Delivery=="Normal" & !is.na(Antibiotics_birth_mother)) 
plot_ordination(antibiotics_mother_filter2,ordination=ordinate(antibiotics_mother_filter2, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(Antibiotics_birth_mother))) +
  stat_ellipse(aes(color=factor(Antibiotics_birth_mother), fill=factor(Antibiotics_birth_mother)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="Antibiotics administered to mother at birth - normal delivery only (P=0.003)     ",breaks=c("1","0"), labels=c("Yes", "No")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.4%)", x="PC1  (23.5%)")
ggsave("plots/antibioticsmother_normalonly_PCoA_subtypes.pdf", width = 11, height =6)

#Antibiotics administered to child, ever
antibiotics_child_1yr_filter <- subset_samples(CRISPRs, !is.na(Antibiotics_ever_1yr_child)) 
plot_ordination(antibiotics_child_1yr_filter,ordination=ordinate(antibiotics_child_1yr_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(Antibiotics_ever_1yr_child))) +
  stat_ellipse(aes(color=factor(Antibiotics_ever_1yr_child), fill=factor(Antibiotics_ever_1yr_child)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="Antibiotics administered to child up until 1yr (P=0.012)    ",breaks=c("1","0"), labels=c("Yes", "No")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=15)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.8%)", x="PC1  (23.1%)")
ggsave("plots/antibioticchildever1yr_PCoA_subtypes.pdf", width = 11, height =6)

#pamcluster
pamcluster_filter <- subset_samples(CRISPRs, !is.na(pamcluster)) 
plot_ordination(pamcluster_filter, ordination=ordinate(pamcluster_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(pamcluster))) +
  stat_ellipse(aes(color=factor(pamcluster), fill=factor(pamcluster)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas system subtypes in children with different types of PAM clusters (P<0.001)      ",breaks=c("2","1"), labels=c("Mature", "Immature")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (18.2%)", x="PC1  (23.3%)")
ggsave("plots/pamcluster_PCoA_subtypes.pdf", width = 11, height =6)

#Asthma 5yr ever
asthma_ever_filter <- subset_samples(CRISPRs, !is.na(j45_5yr_ever)) 
plot_ordination(asthma_ever_filter,ordination=ordinate(asthma_ever_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(j45_5yr_ever))) +
  stat_ellipse(aes(color=factor(j45_5yr_ever), fill=factor(j45_5yr_ever)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas system subtypes in children with asthma diagnosed between age 0-5 (P=0.425)     ",breaks=c("1","0"), labels=c("Yes", "No")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.1%)", x="PC1  (23.3%)")
ggsave("plots/asthma5yrever_PCoA_subtypes.pdf", width = 11, height =6)

#Asthma 5yr cross
asthma_cross_filter <- subset_samples(CRISPRs, !is.na(j45_5yr_cross)) 
plot_ordination(asthma_cross_filter,ordination=ordinate(asthma_cross_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(j45_5yr_cross))) +
  stat_ellipse(aes(color=factor(j45_5yr_cross), fill=factor(j45_5yr_cross)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas system subtypes for children diagnosed with asthma at age 5 (P=0.099)       ",breaks=c("1","0"), labels=c("Yes", "No")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.1%)", x="PC1  (23.3%)")
ggsave("plots/asthma5yrcross_PCoA_subtypes.pdf", width = 11, height =6)

#Asthmatic mother:
asthma_mother_filter <- subset_samples(CRISPRs, !is.na(Asthmatic_mother)) 
plot_ordination(asthma_mother_filter,ordination=ordinate(asthma_mother_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(Asthmatic_mother))) +
  stat_ellipse(aes(color=factor(Asthmatic_mother), fill=factor(Asthmatic_mother)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas system subtypes for children with asthmatic mothers (P=0.218)       ",breaks=c("1","0"), labels=c("Yes", "No")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.7%)", x="PC1  (23.2%)")
ggsave("plots/asthmaticmother_PCoA_subtypes.pdf", width = 11, height =6)

#other analyses:
#Fishoil
oil_filter <- subset_samples(CRISPRs, !is.na(Supplement_oil)) 
plot_ordination(oil_filter,ordination=ordinate(oil_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(Supplement_oil))) +
  stat_ellipse(aes(color=factor(Supplement_oil), fill=factor(Supplement_oil)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas system subtypes in children receiving fishoil or placebo (P=0.651)     ",breaks=c("Fiskeolie","Placebo"), labels=c("Oil", "Placebo")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.7%)", x="PC1  (23.2%)")
ggsave("plots/fishoil_PCoA_subtypes.pdf", width = 11, height =6)

#Delivery
delivery_filter <- subset_samples(CRISPRs, !is.na(Delivery)) 
plot_ordination(delivery_filter,ordination=ordinate(delivery_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(Delivery))) +
  stat_ellipse(aes(color=factor(Delivery), fill=factor(Delivery)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas systems in children born with different methods (P=0.027)       ",breaks=c("Acute sectio","Normal", "Planned sectio"), labels=c("Acute Sectio","Normal", "Planned sectio")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.7%)", x="PC1  (23.2%)")
ggsave("plots/delivery_PCoA_subtypes.pdf", width = 11, height =6)

#Rural/Urban status
ruralurban_filter <- subset_samples(CRISPRs, !is.na(RuralUrbanStatus)) 
plot_ordination(ruralurban_filter,ordination=ordinate(ruralurban_filter, "PCoA", "bray"))+
  theme_bw()+geom_point(aes(color=factor(RuralUrbanStatus))) +
  stat_ellipse(aes(color=factor(RuralUrbanStatus), fill=factor(RuralUrbanStatus)), geom="polygon", level=0.80, size=1, alpha=0.2)+
  scale_color_discrete(name="CRISPR-Cas systems in children with different home placement (P=0.774)       ",breaks=c("Rural","Urban"), labels=c("Rural","Urban")) +
  guides(fill=FALSE) +
  theme(text=element_text(size=13)) +
  theme(legend.position="top") +
  labs(y="PC2  (17.7%)", x="PC1  (23.2%)")
ggsave("plots/ruralurban_PCoA_subtypes.pdf", width = 11, height =6)

