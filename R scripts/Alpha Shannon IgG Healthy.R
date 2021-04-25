## Load packages
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

## in case you need to install the programs
BiocManager::install("ggpubr")



## prune otus that are not present or singletons, From 10928 taxa we went down to 1841 taxa
OTU_HealthyIgG

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 85 samples ]
#sample_data() Sample Data:       [ 85 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]

OTU_HealthyIgGpruned <- prune_taxa(taxa_sums(OTU_HealthyIgG) >= 1, OTU_HealthyIgG)
OTU_HealthyIgGpruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1795 taxa and 85 samples ]
#sample_data() Sample Data:       [ 85 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1795 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1795 tips and 835 internal nodes ]


sort(taxa_sums(OTU_HealthyIgGpruned))

OTU_HealthyIgGpruned_RA1<-rarefy_even_depth(OTU_HealthyIgGpruned, sample.size = 1000, rngseed = 711, replace=FALSE, trimOTUs = TRUE)

#`set.seed(711)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(711); .Random.seed` for the full vector
#...
#11 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
#  sa647sa288sa393sa434sa207
#...
#227OTUs were removed because they are no longer 
#present in any sample after random subsampling

#...

OTU_HealthyIgGpruned_RA1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1568 taxa and 74 samples ]
#sample_data() Sample Data:       [ 74 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1568 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1568 tips and 733 internal nodes ]


## Gives alpha diversity per sample inverse_simpson, gini_simpson, shannon and Fisher

tabOTU_HealthyIgGpruned_RA1 <- alpha(OTU_HealthyIgGpruned_RA1, index = "shannon")
kable(head(tabOTU_HealthyIgGpruned_RA1))

write.table(sample_data (tabOTU_HealthyIgGpruned_RA1), file = "HealthyIgGtShannon.csv", quote = FALSE, sep = "\t", row.names = TRUE)
