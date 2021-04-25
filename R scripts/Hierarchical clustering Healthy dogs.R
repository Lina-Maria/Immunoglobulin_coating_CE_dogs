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

## subset only HealthyIgA


OTU_HealthyIgA <- subset_samples(OTU_Healthy, Subpopulation%in%c("IgA"))
OTU_HealthyIgA

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 87 samples ]
#sample_data() Sample Data:       [ 87 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]

# Correct order for plot
testHealthyIgA <- sample_data(OTU_HealthyIgA)$Population
testHealthyIgA <- factor(testHealthyIgA, levels = c("Total_Bacteria_IgA","Water_presort_IgA","IgA+","IgA-"))
sample_data(OTU_HealthyIgA)$Population <- testHealthyIgA

## subset only HealthyIgA+

OTU_HealthyIgAp <- subset_samples(OTU_HealthyIgA, Population%in%c("IgA+"))
OTU_HealthyIgAp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

# Correct order for plot
testHealthyIgAp <- sample_data(OTU_HealthyIgAp)$Visit
testHealthyIgAp <- factor(testHealthyIgAp, levels = c("1","2"))
sample_data(OTU_HealthyIgAp)$Visit <- testHealthyIgAp


OTU_HealthyIgAppruned <- prune_taxa(taxa_sums(OTU_HealthyIgAp) > 0, OTU_HealthyIgAp)
OTU_HealthyIgAppruned
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 583 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 583 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 583 tips and 382 internal nodes ]


#Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(OTU_HealthyIgAppruned, function(x){x / sum(x)})
phyloseq::otu_table(OTU_HealthyIgAppruned)[1:20, 1:20]

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