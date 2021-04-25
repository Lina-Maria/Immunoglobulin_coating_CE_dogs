library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

## in case you need to install the programs
BiocManager::install("ggpubr")



## prune otus that are not present or singletons, From 10928 taxa we went down to 1841 taxa
OTU_HealthyIgA

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 87 samples ]
#sample_data() Sample Data:       [ 87 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

OTU_HealthyIgApruned <- prune_taxa(taxa_sums(OTU_HealthyIgA) >= 1, OTU_HealthyIgA)
OTU_HealthyIgApruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1687 taxa and 87 samples ]
#sample_data() Sample Data:       [ 87 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 1687 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1687 tips and 885 internal nodes ]

sort(taxa_sums(OTU_HealthyIgApruned))

OTU_HealthyIgApruned_RA1<-rarefy_even_depth(OTU_HealthyIgApruned, sample.size = 1000, rngseed = 711, replace=FALSE, trimOTUs = TRUE)

#`set.seed(711)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(711); .Random.seed` for the full vector
#...
#10 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
#  sa272sa429sa219sa121sa152
#...
#215OTUs were removed because they are no longer 
#present in any sample after random subsampling

OTU_HealthyIgApruned_RA1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1472 taxa and 77 samples ]
#sample_data() Sample Data:       [ 77 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1472 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1472 tips and 784 internal nodes ]

## Gives alpha diversity per sample inverse_simpson, gini_simpson, shannon and Fisher

tabOTU_HealthyIgApruned_RA1 <- alpha(OTU_HealthyIgApruned_RA1, index = "shannon")
kable(head(tabOTU_HealthyIgApruned_RA1))

write.table(sample_data (tabOTU_HealthyIgApruned_RA1), file = "HealthyIgAtShannon.csv", quote = FALSE, sep = "\t", row.names = TRUE)
