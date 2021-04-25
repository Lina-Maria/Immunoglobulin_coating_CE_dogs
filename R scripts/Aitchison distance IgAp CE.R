## Based on Introduction to the Statistical Analysis of Microbiome Data in R by Nicholas Ollberding.

## Load packages
library(microbiome); packageVersion("microbiome")   
#[1] ‘1.10.0’

## CE IgA positive

OTU_CEIgAp <- subset_samples(OTU_CEIgA, Population%in%c("IgA+"))
OTU_CEIgAp

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]


OTU_CEIgAppruned <- prune_taxa(taxa_sums(OTU_CEIgAp) >= 1, OTU_CEIgAp)
OTU_CEIgAppruned

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1140 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1140 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1140 tips and 582 internal nodes ]

(OTU_CEIgAppruned_clr <- microbiome::transform(OTU_CEIgAppruned, "clr")) 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1140 taxa and 46 samples ]
#sample_data() Sample Data:       [ 46 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 1140 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1140 tips and 582 internal nodes ]

phyloseq::otu_table(OTU_CEIgAppruned)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#sa268 sa20 sa137 sa384 sa617
#134ee0689d6529c640e8c4f27f62188e     0    0     0     0     0
#d49b91a991a0bc3f1de4bc7572ffd19d     0    0     0     0     0
#2bbc2d8ec116b9fca2dd50b6515e5461     0    0     0     0     0
#8490a86ce4743c294c1ee4aa4d4701f8     0    0     0     0     0
#8cabd34ed828cc929166c0181ba9a905     0    0     0     0     0

phyloseq::otu_table(OTU_CEIgAppruned_clr)[1:5, 1:5]

#OTU Table:          [5 taxa and 5 samples]
#taxa are rows
#sa268        sa20      sa137      sa384       sa617
#134ee0689d6529c640e8c4f27f62188e -0.1105522 -0.06052073 -0.1429633 -0.1071012 -0.08066729
#d49b91a991a0bc3f1de4bc7572ffd19d -0.1105522 -0.06052073 -0.1429633 -0.1071012 -0.08066729
#2bbc2d8ec116b9fca2dd50b6515e5461 -0.1105522 -0.06052073 -0.1429633 -0.1071012 -0.08066729
#8490a86ce4743c294c1ee4aa4d4701f8 -0.1105522 -0.06052073 -0.1429633 -0.1071012 -0.08066729
#8cabd34ed828cc929166c0181ba9a905 -0.1105522 -0.06052073 -0.1429633 -0.1071012 -0.08066729

OTU_CEIgAppruned_ord <- phyloseq::ordinate(OTU_CEIgAppruned_clr, "RDA")

phyloseq::plot_scree(OTU_CEIgAppruned_ord) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(OTU_CEIgAppruned_ord$CA$eig)  

#PC1      PC2      PC3      PC4      PC5      PC6 
#50.03794 33.26226 30.60183 26.36437 24.92560 24.28675 

sapply(OTU_CEIgAppruned_ord$CA$eig[1:5], function(x) x / sum(OTU_CEIgAppruned_ord$CA$eig))

#PC1        PC2        PC3        PC4        PC5 
#0.07083921 0.04708971 0.04332331 0.03732430 0.03528742 

clr1 <- OTU_CEIgAppruned_ord$CA$eig[1] / sum(OTU_CEIgAppruned_ord$CA$eig)
clr2 <- OTU_CEIgAppruned_ord$CA$eig[2] / sum(OTU_CEIgAppruned_ord$CA$eig)
phyloseq::plot_ordination(OTU_CEIgAppruned, OTU_CEIgAppruned_ord, type="samples", color="Stage") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Stage), linetype = 2)


clr_dist_matrix <- phyloseq::distance(OTU_CEIgAppruned_clr, method = "euclidean") 
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(OTU_CEIgAppruned_clr)$Stage)

#Call:
#  vegan::adonis(formula = clr_dist_matrix ~ phyloseq::sample_data(OTU_CEIgAppruned_clr)$Stage) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#phyloseq::sample_data(OTU_CEIgAppruned_clr)$Stage  1       595  594.75 0.83898 0.01871  0.974
#Residuals                                         44     31191  708.90         0.98129       
#Total                                             45     31786                 1.00000


dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(OTU_CEIgAppruned_clr)$Stage)
dispr

#Homogeneity of multivariate dispersions

#Call: vegan::betadisper(d = clr_dist_matrix, group =
#                          phyloseq::sample_data(OTU_CEIgAppruned_clr)$Stage)

#No. of Positive Eigenvalues: 45
#No. of Negative Eigenvalues: 0

#Average distance to median:
#  Active Remission 
#26.20     24.97 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 45 eigenvalues)
#PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
#2252  1497  1377  1186  1122  1093  1051  1031 


plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

permutest(dispr)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     1  16.79  16.788 0.9438    999  0.345
#Residuals 44 782.68  17.788   
