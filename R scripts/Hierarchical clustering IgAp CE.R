library(vegan); packageVersion("vegan")     
#[1] ‘2.5.6’

## subset only CE

OTU_CE <- subset_samples(OTU_FaecesWater, Group%in%c("ChronicEnteropathy"))
OTU_CE

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 347 samples ]
#sample_data() Sample Data:       [ 347 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]

# Correct order for plot
testCE <- sample_data(OTU_CE)$Population
testCE <- factor(testCE, levels = c("Total_Bacteria_IgA","Water_presort_IgA","IgA+","IgA-","Total_Bacteria_IgG","Water_presort_IgG","IgG+","IgG-"))
sample_data(OTU_CE)$Population <- testCE


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

OTU_CEIgAp <- subset_samples(OTU_CEIgA, Population%in%c("IgA+"))
OTU_CEIgAp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

OTU_CEIgAppruned <- prune_taxa(taxa_sums(OTU_CEIgAp) >= 1, OTU_CEIgAp)
OTU_CEIgAppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1140 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 1140 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1140 tips and 582 internal nodes ]


#Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(OTU_CEIgAppruned, function(x){x / sum(x)})
phyloseq::otu_table(OTU_CEIgAppruned)[1:20, 1:20]

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
colorCode <- c(Noodle = "red", Winnie_Lee = "blue", Mouse = "blueviolet", Toby_Walker = "brown4", Toby_Hall = "darkblue", Millie_Abbott = "darkgoldenrod4", Cruze = "darkgreen", Jetta_Krop = "darkmagenta", Millie_Hurley = "darkorange1", TimTam = "darkorchid4", Boof_Hedges = "darkred", Pissqueek = "deeppink3", Bonnie = "deepskyblue4", Mako_PDU = "firebrick1", Frankie = "green1", Shiraz_McPherson = "black", Oxley = "yellow1", Bridgette = "burlywood3", Bentley_Smith = "cadetblue1", Buckley = "cornflowerblue")

labels_colors(ward) <- colorCode[meta$Name][order.dendrogram(ward)]



#Plot
plot(ward)