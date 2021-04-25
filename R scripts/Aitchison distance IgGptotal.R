## Based on Introduction to the Statistical Analysis of Microbiome Data in R by Nicholas Ollberding.

## Load packages
library(microbiome); packageVersion("microbiome")   
#[1] ‘1.10.0’

# subset only TotalIgG+

OTU_TotalIgGp <- subset_samples(OTU_FaecesWater, Population%in%c("IgG+"))
OTU_TotalIgGp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 67 samples ]
#sample_data() Sample Data:       [ 67 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]

OTU_TotalIgGppruned <- prune_taxa(taxa_sums(OTU_TotalIgGp) > 0, OTU_TotalIgGp)
OTU_TotalIgGppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1596 taxa and 67 samples ]
#sample_data() Sample Data:       [ 67 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1596 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1596 tips and 734 internal nodes ]

(OTU_TotalIgGppruned_clr <- microbiome::transform(OTU_TotalIgGppruned, "clr")) 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1596 taxa and 67 samples ]
#sample_data() Sample Data:       [ 67 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1596 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1596 tips and 734 internal nodes ]

phyloseq::otu_table(OTU_TotalIgGppruned)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#sa122 sa432 sa139 sa336 sa45
#ac980b106589bff32bbd91f5196a1251     0     0     0     0    0
#e73017fd7adaa63ab1294170b0cb5264     0     0     0     0    0
#c325db87106491b3dcd98b0ace3a4e2d     0     0     0     0    0
#feec183d85ed9f5db70aa84ca2c649d5     0     0     0     0    0
#239436d77092b03ab44e7f2a336159e4     0     0     0     0    0

phyloseq::otu_table(OTU_TotalIgGppruned_clr)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#sa122      sa432       sa139      sa336      sa45
#ac980b106589bff32bbd91f5196a1251 -0.1124155 -0.1531688 -0.06088343 -0.0633759 -0.164088
#e73017fd7adaa63ab1294170b0cb5264 -0.1124155 -0.1531688 -0.06088343 -0.0633759 -0.164088
#c325db87106491b3dcd98b0ace3a4e2d -0.1124155 -0.1531688 -0.06088343 -0.0633759 -0.164088
#feec183d85ed9f5db70aa84ca2c649d5 -0.1124155 -0.1531688 -0.06088343 -0.0633759 -0.164088
#239436d77092b03ab44e7f2a336159e4 -0.1124155 -0.1531688 -0.06088343 -0.0633759 -0.164088

OTU_TotalIgGppruned_ord <- phyloseq::ordinate(OTU_TotalIgGppruned_clr, "RDA")

phyloseq::plot_scree(OTU_TotalIgGppruned_ord) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(OTU_TotalIgGppruned_ord$CA$eig)  

# PC1      PC2      PC3      PC4      PC5      PC6 
#50.71695 34.82163 32.87103 29.06417 26.75743 24.23036 

sapply(OTU_TotalIgGppruned_ord$CA$eig[1:5], function(x) x / sum(OTU_TotalIgGppruned_ord$CA$eig))

#PC1        PC2        PC3        PC4        PC5 
#0.05426222 0.03725577 0.03516882 0.03109585 0.02862785

clr1 <- OTU_TotalIgGppruned_ord$CA$eig[1] / sum(OTU_TotalIgGppruned_ord$CA$eig)
clr2 <- OTU_TotalIgGppruned_ord$CA$eig[2] / sum(OTU_TotalIgGppruned_ord$CA$eig)
phyloseq::plot_ordination(OTU_TotalIgGppruned, OTU_TotalIgGppruned_ord, type="samples", color="Group") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Group), linetype = 2)


clr_dist_matrix <- phyloseq::distance(OTU_TotalIgGppruned_clr, method = "euclidean") 
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(OTU_TotalIgGppruned_clr)$Group)

#Call:
#  vegan::adonis(formula = clr_dist_matrix ~ phyloseq::sample_data(OTU_TotalIgGppruned_clr)$Group) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#phyloseq::sample_data(OTU_TotalIgGppruned_clr)$Group  1       961  960.90  1.0285 0.01558  0.322
#Residuals                                            65     60727  934.26         0.98442       
#Total                                                66     61688                 1.00000       

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(OTU_TotalIgGppruned_clr)$Group)
dispr

#Homogeneity of multivariate dispersions

#Call: vegan::betadisper(d = clr_dist_matrix, group =
#                          phyloseq::sample_data(OTU_TotalIgGppruned_clr)$Group)

#No. of Positive Eigenvalues: 66
#No. of Negative Eigenvalues: 0

#Average distance to median:
#  ChronicEnteropathy            Healthy 
#30.10              28.62 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 66 eigenvalues)
#PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
#3347  2298  2169  1918  1766  1599  1542  1527 

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

permutest(dispr)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     1   31.56  31.555 1.0776        
#Residuals 65 1903.41  29.283

