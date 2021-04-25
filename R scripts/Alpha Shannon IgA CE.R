library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

## in case you need to install the programs
BiocManager::install("ggpubr")



## prune otus that are not present or singletons, From 10928 taxa we went down to 1841 taxa
## subset only CEIgA

OTU_CEIgA <- subset_samples(OTU_CE, Subpopulation%in%c("IgA"))
OTU_CEIgA

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 175 samples ]
#sample_data() Sample Data:       [ 175 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]

# Correct order for plot
testCEIgA <- sample_data(OTU_CEIgA)$Population
testCEIgA <- factor(testCEIgA, levels = c("Total_Bacteria_IgA","Water_presort_IgA","IgA+","IgA-"))
sample_data(OTU_CEIgA)$Population <- testCEIgA


OTU_CEIgApruned <- prune_taxa(taxa_sums(OTU_CEIgA) >= 1, OTU_CEIgA)
OTU_CEIgApruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3450 taxa and 175 samples ]
#sample_data() Sample Data:       [ 175 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 3450 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 3450 tips and 1441 internal nodes ]


sort(taxa_sums(OTU_CEIgApruned))

OTU_CEIgApruned_RA1<-rarefy_even_depth(OTU_CEIgApruned, sample.size = 1000, rngseed = 711, replace=FALSE, trimOTUs = TRUE)

#`set.seed(711)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(711); .Random.seed` for the full vector
#...
#25 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
#  sa335sa74sa403sa307sa75
#...
#431OTUs were removed because they are no longer 
#present in any sample after random subsampling

#...
OTU_CEIgApruned_RA1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3019 taxa and 150 samples ]
#sample_data() Sample Data:       [ 150 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 3019 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 3019 tips and 1278 internal nodes ]

## Gives alpha diversity per sample inverse_simpson, gini_simpson, shannon and Fisher

tabOTU_CEIgApruned_RA1 <- alpha(OTU_CEIgApruned_RA1, index = "shannon")
kable(head(tabOTU_CEIgApruned_RA1))

write.table(sample_data (tabOTU_CEIgApruned_RA1), file = "IgACEShannon.csv", quote = FALSE, sep = "\t", row.names = TRUE)