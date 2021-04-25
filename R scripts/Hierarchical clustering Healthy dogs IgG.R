library(vegan); packageVersion("vegan")     
#[1] ‘2.5.6’

## subset only Healthy

OTU_Healthy <- subset_samples(OTU_FaecesWater, Group%in%c("Healthy"))
OTU_Healthy

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 172 samples ]
#sample_data() Sample Data:       [ 172 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

## subset only HealthyIgG

OTU_HealthyIgG <- subset_samples(OTU_Healthy, Subpopulation%in%c("IgG"))
OTU_HealthyIgG

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 85 samples ]
#sample_data() Sample Data:       [ 85 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

# Correct order for plot
testHealthyIgG <- sample_data(OTU_HealthyIgG)$Population
testHealthyIgG <- factor(testHealthyIgG, levels = c("Total_Bacteria_IgG","Water_presort_IgG","IgG+","IgG-"))
sample_data(OTU_HealthyIgG)$Population <- testHealthyIgG

## subset only HealthyIgG+

OTU_HealthyIgGp <- subset_samples(OTU_HealthyIgG, Population%in%c("IgG+"))
OTU_HealthyIgGp



#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 21 samples ]
#sample_data() Sample Data:       [ 21 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]

# Correct order for plot
testHealthyIgGp <- sample_data(OTU_HealthyIgGp)$Visit
testHealthyIgGp <- factor(testHealthyIgGp, levels = c("1","2"))
sample_data(OTU_HealthyIgGp)$Visit <- testHealthyIgGp


OTU_HealthyIgGppruned <- prune_taxa(taxa_sums(OTU_HealthyIgGp) > 0, OTU_HealthyIgGp)
OTU_HealthyIgGppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 570 taxa and 21 samples ]
#sample_data() Sample Data:       [ 21 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 570 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 570 tips and 343 internal nodes ]


#Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(OTU_HealthyIgGppruned, function(x){x / sum(x)})
phyloseq::otu_table(OTU_HealthyIgGppruned)[1:20, 1:20]

phyloseq::otu_table(ps_rel_abund)[1:20, 1:20]

#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]

#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "average"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(Jedd = "red", Savannah_Gillies = "blue", Fergus =  "blueviolet", Polly = "brown4", Ellie = "darkblue", Finnegan = "darkgoldenrod4", Diesel = "darkgreen",  Sophie = "darkmagenta", Munstar = "darkorange1", Skyler = "darkorchid4", Electre = "darkred")

labels_colors(ward) <- colorCode[meta$Name][order.dendrogram(ward)]


#Plot
plot(ward)