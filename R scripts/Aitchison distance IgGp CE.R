## Based on Introduction to the Statistical Analysis of Microbiome Data in R by Nicholas Ollberding.

## Load packages
library(microbiome); packageVersion("microbiome")   
#[1] ‘1.10.0’


## CE IgG positive

OTU_CEIgGp <- subset_samples(OTU_CEIgG, Population%in%c("IgG+"))
OTU_CEIgGp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]


OTU_CEIgGppruned <- prune_taxa(taxa_sums(OTU_CEIgGp) >= 1, OTU_CEIgGp)
OTU_CEIgGppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1163 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1163 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1163 tips and 574 internal nodes ]

(OTU_CEIgGppruned_clr <- microbiome::transform(OTU_CEIgGppruned, "clr")) 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1163 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1163 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1163 tips and 574 internal nodes ]

phyloseq::otu_table(OTU_CEIgGppruned)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#s                                 sa139 sa336 sa45 sa636 sa76
#ac980b106589bff32bbd91f5196a1251     0     0    0     0    0
#c325db87106491b3dcd98b0ace3a4e2d     0     0    0     0    0
#feec183d85ed9f5db70aa84ca2c649d5     0     0    0     0    0
#239436d77092b03ab44e7f2a336159e4     0     0    0     0    0
#4098ce3a8a9cb92b60cdb74fe4c4ee25     0     0    0     0    0

phyloseq::otu_table(OTU_CEIgGppruned_clr)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#                                       sa139       sa336       sa45      sa636        sa76
#ac980b106589bff32bbd91f5196a1251 -0.07953105 -0.08254767 -0.2155238 -0.1265445 -0.08921542
#c325db87106491b3dcd98b0ace3a4e2d -0.07953105 -0.08254767 -0.2155238 -0.1265445 -0.08921542
#feec183d85ed9f5db70aa84ca2c649d5 -0.07953105 -0.08254767 -0.2155238 -0.1265445 -0.08921542
#239436d77092b03ab44e7f2a336159e4 -0.07953105 -0.08254767 -0.2155238 -0.1265445 -0.08921542
#4098ce3a8a9cb92b60cdb74fe4c4ee25 -0.07953105 -0.08254767 -0.2155238 -0.1265445 -0.08921542

OTU_CEIgGppruned_ord <- phyloseq::ordinate(OTU_CEIgGppruned_clr, "RDA")

phyloseq::plot_scree(OTU_CEIgGppruned_ord) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(OTU_CEIgGppruned_ord$CA$eig)  

#PC1      PC2      PC3      PC4      PC5      PC6 
#64.95347 36.43190 35.56720 33.77562 30.99153 30.08641 

sapply(OTU_CEIgGppruned_ord$CA$eig[1:5], function(x) x / sum(OTU_CEIgGppruned_ord$CA$eig))

#PC1        PC2        PC3        PC4        PC5 
#0.07083921 0.04708971 0.04332331 0.03732430 0.03528742 

clr1 <- OTU_CEIgGppruned_ord$CA$eig[1] / sum(OTU_CEIgGppruned_ord$CA$eig)
clr2 <- OTU_CEIgGppruned_ord$CA$eig[2] / sum(OTU_CEIgGppruned_ord$CA$eig)
phyloseq::plot_ordination(OTU_CEIgGppruned, OTU_CEIgGppruned_ord, type="samples", color="Stage") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Stage), linetype = 2)


clr_dist_matrix <- phyloseq::distance(OTU_CEIgGppruned_clr, method = "euclidean") 
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(OTU_CEIgGppruned_clr)$Stage)

#CCall:
#vegan::adonis(formula = clr_dist_matrix ~ phyloseq::sample_data(OTU_CEIgGppruned_clr)$Stage) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#phyloseq::sample_data(OTU_CEIgGppruned_clr)$Stage  1       702  701.77 0.79876 0.01783  0.998
#Residuals                                         44     38657  878.57         0.98217       
#Total                                             45     39359                 1.00000    


dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(OTU_CEIgGppruned_clr)$Stage)
dispr

#Homogeneity of multivariate dispersions

#Call: vegan::betadisper(d = clr_dist_matrix, group =
#                          phyloseq::sample_data(OTU_CEIgGppruned_clr)$Stage)

#No. of Positive Eigenvalues: 45
#No. of Negative Eigenvalues: 0

#Average distance to median:
#  Active Remission 
#29.35     27.43 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 45 eigenvalues)
#PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
#2923  1639  1601  1520  1395  1354  1308  1301 


plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

permutest(dispr)

#permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     1   40.72  40.717 1.7214    999  0.234
#Residuals 44 1040.71  23.653        