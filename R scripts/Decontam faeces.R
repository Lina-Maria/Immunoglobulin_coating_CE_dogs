install.packages("devtools")
library(devtools)
devtools::install_github("benjjneb/decontam")

library(decontam)
packageVersion("decontam")
?isContaminant

## pHYSEQ1 IN GENENERAL
#Inspect Library Sizes
#Let’s take a quick first look at the library sizes (i.e. the number of reads) in each sample, as a function of whether that sample was a true positive sample or a negative control:

head(sample_data(Physeq1))
df <- as.data.frame(sample_data(Physeq1))
df$LibrarySize <- sample_sums(Physeq1)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


#Identify Contaminants - Prevalence
#In this method, the prevalence (presence/absence across samples) of each sequence feature in true positive samples is compared to the prevalence in negative controls to identify contaminants.
#In our phyloseq object, "Sample_or_Control" is the sample variable that holds the negative control information. We’ll summarize that data as a logical variable, with TRUE for control samples, as that is the form required by isContaminant.
## pHYSEQ1 TAKING INTO ACCOUNT ONLY THE FAECES POSITIVE CONTROL AND RELEVANT NEGATIVE CONTROLS FOR THESE SAMPLES AND THE BATCHES: IN THIS CASE RUN: THRESHOLD 0.1
sample_data(Physeq1)$is.neg <- sample_data(Physeq1)$Sample_or_Control == "Control Sample"
contamdf1.run <- isContaminant(Physeq1, method="prevalence", neg="is.neg", batch = "Run")
table(contamdf1.run$contaminant)

#FALSE  TRUE 
#10830    98 

head(which(contamdf1.run$contaminant))
#[1]  212  898  904 1074 1315 1362

write.table(contamdf1.run, file = "contamdf1.run.csv", quote = FALSE, sep = "\t", row.names = TRUE)

print(which(contamdf1.run$contaminant)) ### THESE IS THE LIST OF ALL OTUS THAT ARE CONTAMINANT THE POSITION NUMBER IN THE TABLE

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(Physeq1, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf1.run$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

## pHYSEQ1 TAKING INTO ACCOUNT ONLY THE FAECES POSITIVE CONTROL AND RELEVANT NEGATIVE CONTROLS FOR THESE SAMPLES AND THE BATCHES: IN THIS CASE RUN: THRESHOLD 0.5

contamdf1.prev05.run <- isContaminant(Physeq1, method="prevalence", neg="is.neg", threshold=0.5, batch = "Run")
table(contamdf1.prev05.run$contaminant)

#FALSE  TRUE 
#10609   319

head(which(contamdf1.prev05.run$contaminant))
#[1]  37  99 142 181 212 241

write.table(contamdf1.prev05.run, file = "contamdf1.prev05.run.csv", quote = FALSE, sep = "\t", row.names = TRUE)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(Physeq1, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf1.prev05.run$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


## How to do the histogram
hist(contamdf1.prev05.run$p, 100, ylim = c(0,300), xlim = c(0,1))
hist(contamdf1.run$p, 100, ylim = c(0,300), xlim = c(0,1))

Physeq1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10928 taxa and 544 samples ]
#sample_data() Sample Data:       [ 544 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10928 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10928 tips and 3566 internal nodes ]


#Now that we have identified likely contaminants, let’s remove them from the phyloseq object: THRESHOLD 0.5

Physeq1.noncontam1 <- prune_taxa(!contamdf1.prev05.run$contaminant, Physeq1)
Physeq1.noncontam1

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10609 taxa and 544 samples ]
#sample_data() Sample Data:       [ 544 samples by 30 sample variables ]
#tax_table()   Taxonomy Table:    [ 10609 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 10609 tips and 3547 internal nodes ]