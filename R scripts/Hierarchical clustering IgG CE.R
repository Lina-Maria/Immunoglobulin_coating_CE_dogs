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

## subset only CEIgG+

OTU_CEIgGp <- subset_samples(OTU_CEIgG, Population%in%c("IgG+"))
OTU_CEIgGp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]


OTU_CEIgGppruned <- prune_taxa(taxa_sums(OTU_CEIgGp) >= 1, OTU_CEIgGp)
OTU_CEIgGppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1163 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1163 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1163 tips and 574 internal nodes ]


#Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(OTU_CEIgGppruned, function(x){x / sum(x)})
phyloseq::otu_table(OTU_CEIgGppruned)[1:20, 1:20]

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