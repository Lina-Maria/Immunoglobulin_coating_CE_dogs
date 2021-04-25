library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

## in case you need to install the programs
BiocManager::install("ggpubr")



## prune otus that are not present or singletons, From 10928 taxa we went down to 1841 taxa
## subset only CEIgG

OTU_CEIgG <- subset_samples(OTU_CE, Subpopulation%in%c("IgG"))
OTU_CEIgG

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 172 samples ]
#sample_data() Sample Data:       [ 172 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

# Correct order for plot
testCEIgG <- sample_data(OTU_CEIgG)$Population
testCEIgG <- factor(testCEIgG, levels = c("Total_Bacteria_IgG","Water_presort_IgG","IgG+","IgG-"))
sample_data(OTU_CEIgG)$Population <- testCEIgG



OTU_CEIgGpruned <- prune_taxa(taxa_sums(OTU_CEIgG) >= 1, OTU_CEIgG)
OTU_CEIgGpruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3604 taxa and 172 samples ]
#sample_data() Sample Data:       [ 172 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 3604 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 3604 tips and 1427 internal nodes ]

sort(taxa_sums(OTU_CEIgGpruned))

OTU_CEIgGpruned_RA1<-rarefy_even_depth(OTU_CEIgGpruned, sample.size = 1000, rngseed = 711, replace=FALSE, trimOTUs = TRUE)


#`set.seed(711)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(711); .Random.seed` for the full vector
#...
#16 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
#  sa64sa421sa296sa490sa17
#...
#401OTUs were removed because they are no longer 
#present in any sample after random subsampling

OTU_CEIgGpruned_RA1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3203 taxa and 156 samples ]
#sample_data() Sample Data:       [ 156 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 3203 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 3203 tips and 1294 internal nodes ]

## Gives alpha diversity per sample inverse_simpson, gini_simpson, shannon and Fisher

tabOTU_CEIgGpruned_RA1 <- alpha(OTU_CEIgGpruned_RA1, index = "shannon")
kable(head(tabOTU_CEIgGpruned_RA1))

write.table(sample_data (tabOTU_CEIgGpruned_RA1), file = "IgGCEShannon.csv", quote = FALSE, sep = "\t", row.names = TRUE)
