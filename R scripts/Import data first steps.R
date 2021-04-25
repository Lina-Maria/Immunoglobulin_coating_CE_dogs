## Load packages
library("phyloseq")
library("ggplot2")


otu_table <- read.csv("otu_matrix.csv", sep = ",", row.names = 1)
otu_table <- as.matrix(otu_table)

OTU <- otu_table(otu_table, taxa_are_rows = TRUE)


taxonomy = read.csv("taxonomy1.csv", sep=",", row.names=1)
taxonomy = as.matrix (taxonomy)

TAX = tax_table(taxonomy) 

## the file must be saved in open office calculator and then copy and paste in excel and save as csv. Be careful that sample name contains sa plus number sa125
metadata	=	read.csv("MetadataFCT_CE_merged3.csv", sep=",", row.names=1)
META = sample_data(metadata)


phy_tree = read_tree("tree.nwk")


taxa_names(TAX) 
taxa_names(OTU) 
taxa_names(phy_tree)

sample_names(OTU) 
sample_names(META)

physeq = phyloseq(OTU, TAX, META, phy_tree)
physeq

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10928 taxa and 628 samples ]
#sample_data() Sample Data:       [ 628 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10928 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10928 tips and 3566 internal nodes ]

## if you want to check every component in detail:
##head: only the first rows of information
## tail: only the last rows
##print: everything

## for sample data use:
head(sample_data(physeq))
tail(sample_data(physeq))
print(sample_data(physeq))


write.table(sample_data (physeq), file = "physeqSample.csv", quote = FALSE, sep = "\t", row.names = TRUE)

## for OTU table, replace sample_data for Otu_table and for taxonomy table, replace sample_date for tax_table:
head(otu_table(physeq))
head(tax_table(physeq))


