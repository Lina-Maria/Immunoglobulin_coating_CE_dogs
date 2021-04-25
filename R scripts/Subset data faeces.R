## subset only faeces, water presort and negatives

Physeq1 <- subset_samples(physeq, SampleType%in%c("Water_presort","Faeces","Water_negative","Tris_Buffer","Mock_community","NA","Interun_control","Water_WEHI"))
Physeq1

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10928 taxa and 544 samples ]
#sample_data() Sample Data:       [ 544 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10928 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10928 tips and 3566 internal nodes ]


# subset only faeces and water presort

OTU_FaecesWater <- subset_samples(Physeq1.noncontam1, SampleType%in%c("Water_presort","Faeces"))
OTU_FaecesWater

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 519 samples ]
#sample_data() Sample Data:       [ 519 samples by 31 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]