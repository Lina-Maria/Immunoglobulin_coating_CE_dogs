## Based on Introduction to the Statistical Analysis of Microbiome Data in R by Nicholas Ollberding.

## Load packages
library(microbiome); packageVersion("microbiome")   
#[1] ‘1.10.0’

# subset only TotalIgA+

OTU_TotalIgAp <- subset_samples(OTU_FaecesWater, Population%in%c("IgA+"))
OTU_TotalIgAp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3566 internal nodes ]

OTU_TotalIgAppruned <- prune_taxa(taxa_sums(OTU_TotalIgAp) > 0, OTU_TotalIgAp)
OTU_TotalIgAppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1555 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1555 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1555 tips and 776 internal nodes ]

(OTU_TotalIgAppruned_clr <- microbiome::transform(OTU_TotalIgAppruned, "clr")) 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1555 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1555 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1555 tips and 776 internal nodes ]

phyloseq::otu_table(OTU_TotalIgAppruned)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#sa631 sa14 sa268 sa20 sa137
#4bdbb19e6a410287627ddf2c2166a534     0    0     0    0     0
#d0ad8b8b41bc46299a3bfb67aa3fa8a0     0    0     0    0     0
#134ee0689d6529c640e8c4f27f62188e     0    0     0    0     0
#7305e4cdc630c73e73a9f068b5dd367c     0    0     0    0     0
#332faf523622bf746eb8e243f3e2f1af     0    0     0    0     0

phyloseq::otu_table(OTU_TotalIgAppruned_clr)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#sa631       sa14       sa268        sa20      sa137
#4bdbb19e6a410287627ddf2c2166a534 -0.1227948 -0.1093876 -0.08359953 -0.04590017 -0.1079212
#d0ad8b8b41bc46299a3bfb67aa3fa8a0 -0.1227948 -0.1093876 -0.08359953 -0.04590017 -0.1079212
#134ee0689d6529c640e8c4f27f62188e -0.1227948 -0.1093876 -0.08359953 -0.04590017 -0.1079212
#7305e4cdc630c73e73a9f068b5dd367c -0.1227948 -0.1093876 -0.08359953 -0.04590017 -0.1079212
#332faf523622bf746eb8e243f3e2f1af -0.1227948 -0.1093876 -0.08359953 -0.04590017 -0.1079212

OTU_TotalIgAppruned_ord <- phyloseq::ordinate(OTU_TotalIgAppruned_clr, "RDA")

phyloseq::plot_scree(OTU_TotalIgAppruned_ord) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(OTU_TotalIgAppruned_ord$CA$eig)  

# PC1      PC2      PC3      PC4      PC5      PC6 
#38.18650 28.96693 26.36318 24.86377 22.43039 19.40285

sapply(OTU_TotalIgAppruned_ord$CA$eig[1:5], function(x) x / sum(OTU_TotalIgAppruned_ord$CA$eig))

#PC1        PC2        PC3        PC4        PC5 
#0.05236801 0.03972453 0.03615380 0.03409755 0.03076047

clr1 <- OTU_TotalIgAppruned_ord$CA$eig[1] / sum(OTU_TotalIgAppruned_ord$CA$eig)
clr2 <- OTU_TotalIgAppruned_ord$CA$eig[2] / sum(OTU_TotalIgAppruned_ord$CA$eig)
phyloseq::plot_ordination(OTU_TotalIgAppruned, OTU_TotalIgAppruned_ord, type="samples", color="Group") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Group), linetype = 2)


clr_dist_matrix <- phyloseq::distance(OTU_TotalIgAppruned_clr, method = "euclidean") 
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(OTU_TotalIgAppruned_clr)$Group)

#Call:
#  vegan::adonis(formula = clr_dist_matrix ~ phyloseq::sample_data(OTU_TotalIgAppruned_clr)$Group) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#phyloseq::sample_data(OTU_TotalIgAppruned_clr)$Group  1       724  724.21 0.99305 0.01482  0.517
#Residuals                                            66     48132  729.27         0.98518       
#Total                                                67     48856                 1.00000        

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(OTU_TotalIgAppruned_clr)$Group)
dispr

#Homogeneity of multivariate dispersions

#Call: vegan::betadisper(d = clr_dist_matrix, group =
#                          phyloseq::sample_data(OTU_TotalIgAppruned_clr)$Group)

#No. of Positive Eigenvalues: 67
#No. of Negative Eigenvalues: 0

#Average distance to median:
#  ChronicEnteropathy            Healthy 
#26.68              25.29 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 67 eigenvalues)
#PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
#2558  1941  1766  1666  1503  1300  1267  1187 

#plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
#boxplot(dispr, main = "", xlab = "")

permutest(dispr)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     1   28.59  28.587 1.3767    999  0.249
#Residuals 66 1370.46  20.765     
