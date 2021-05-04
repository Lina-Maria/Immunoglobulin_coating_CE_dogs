# AP. Woodward and L. Martinez-Lopez, 2020-21, University of Melbourne.
# Implements a series of analyses of alpha-diversity, flow-cytometry antibody binding, and antibody-dependent taxon relative abundance.
# Post-processing and plot generation is also implemented.
# At each stage the analytical datasets are exported for reproducibility.


# Conducts analysis of alpha-diversity (Shannon index) from the microbiome dataset including antibody status, disease classification, and stage
#     as predictors.
# Then plots and summarizes the results.

library(lme4)
library(emmeans)
library(MuMIn)
library(ggplot2)
library(brms)
library(performance)
library(tidybayes)
library(ggpubr)


# Conduct a linear mixed model for the log(Shannon) with random intercept for Dog.
#     log-transformation is natural for Shannon index as it is positive-nonzero.
#     The transformation is readily justified by the diagnostic plots.
# The full set of possible interacting terms are used here.

Shannon_mod <- lmer(log(Shannon)~1+Ig*Classification_CE*Stage+(1|Subject), data = antibody_Shannon_deID)

summary(Shannon_mod)
confint(Shannon_mod)

hist(residuals(Shannon_mod))
plot(fitted(Shannon_mod), residuals(Shannon_mod))

r.squaredGLMM(Shannon_mod)

write.csv(antibody_Shannon_deID, file = 'alpha_antibody_data_AW_deID.csv', row.names = FALSE)

# Using the 'emmeans' package, develop predicted expected value of Shannon under each condition, and then generate an interval plot of the results.
#     Also print the 'emmeans' table.

Shannon_post_means <- emmeans(Shannon_mod, specs = c('Ig','Classification_CE','Stage'), type = 'response')

Shannon_post_data <- as.data.frame(Shannon_post_means)
Shannon_post_data$Classification_CE <- factor(Shannon_post_data$Classification_CE, levels = c('Healthy','ARE','DRE','IRE'))
Shannon_post_data$Stage <- factor(Shannon_post_data$Stage, levels = c('Before','After'))
Shannon_post_data$Ig <- factor(Shannon_post_data$Ig, levels = c('Input','IgAp','IgGp'))

shannon_index_plot <- ggplot(data = Shannon_post_data, aes(x = Classification_CE, color = Stage, y = response, ymin = lower.CL, ymax = upper.CL)) + geom_point(shape = 2, position = position_dodge(width = 0.7)) + geom_errorbar(position = position_dodge(width = 0.7)) + facet_wrap(~Ig, nrow = 1) + theme(axis.title.x=element_blank()) + ylab('Estimated Shannon Index') + ylim(c(1.5,4.5))
ggsave(filename = 'shannon_index_plot.tif', plot = shannon_index_plot, device = 'tiff', dpi = 100, width = 7, height = 3)

print(Shannon_post_means)


# Implements binomial hierarchical analyses of flow cytometry data under IgA or IgG binding conditions.
# Models are estimated using the binomial() family in package 'brms' (10.18637/jss.v080.i01).
# Once the models are completed, visualization of posterior proportion positive is generated, and a table of the values is returned.

antibody_data_revised_APW$Classification_CE <- factor(antibody_data_revised_APW$Classification_CE, levels = c('Healthy','ARE','FRE','IRE'))
antibody_data_revised_APW$Stage <- factor(antibody_data_revised_APW$Stage, levels = c('Before','After'))

antibody_dcount_mod1 <- brm(bf(Positive|trials(Total)~1+Antibody*Stage*Classification_CE+(1+Antibody|Subject)), family = binomial(), data = antibody_data_revised_APW_deID, cores = 4, control = list(max_treedepth = 15), iter = 3000)
r2_bayes(antibody_dcount_mod1)

new_antibody_data <- data.frame(Antibody = c('A','A','A','A','A','A','A','A','G','G','G','G','G','G','G','G'), Stage = c('Before','Before','Before','Before','After','After','After','After','Before','Before','Before','Before','After','After','After','After'), Classification_CE = c('FRE','ARE','IRE','Healthy','FRE','ARE','IRE','Healthy','FRE','ARE','IRE','Healthy','FRE','ARE','IRE','Healthy'), Total = 1)

new_antibody_estimates <- fitted(antibody_dcount_mod1, newdata = new_antibody_data, re_formula = NA, probs = c(0.05,0.95))
new_antibody_totals <- cbind(new_antibody_data,new_antibody_estimates)

new_antibody_totals$Stage <- factor(new_antibody_totals$Stage, levels = c('Before','After'))
new_antibody_totals$Classification_CE <- factor(new_antibody_totals$Classification_CE, levels = c('Healthy','ARE','FRE','IRE'))

Antibody_names <- c('Immunoglobulin-A','Immunoglobulin-G')
names(Antibody_names) <- c('A','G')
ggplot(data = new_antibody_totals, aes(x = Stage, y = Estimate, color = Classification_CE)) + geom_point(shape = 2, size = 3, position = (position_dodge(width = 0.8))) + facet_wrap(~Antibody, labeller = labeller(Antibody = Antibody_names)) + geom_errorbar(aes(ymin = (Q5), ymax = (Q95)), position = position_dodge(width = 0.8)) + geom_point(data = antibody_data_revised_APW, size = 1.5, aes(x = Stage, y = Percentage, color = Classification_CE), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.05)) + ylim(c(0,1)) + theme(axis.title.x=element_blank()) + ylab('Posterior Proportion Positive')

print(new_antibody_totals)

write.csv(antibody_data_revised_APW_deID, file = 'flow_cytometry_data_AW_deID.csv', row.names = FALSE)

# Implements multinomial hierarchical analyses of OTU data under IgA or IgG binding conditions.
# Models are estimated using the multinomial() family in package 'brms' (10.18637/jss.v080.i01).
#     Note that this analysis is computationally expensive (about 8 hours runtime on a desktop PC with 16gb RAM and 3.9GHz CPU).
# Once the models are completed, visualization of posterior proportional abundances are generated.
# To examine the effects of IgA and IgG binding, we generate the 'Palm' index (), from the posterior proportional abundances
#     rather than the observed OTU. We determine the Ig+:Ig- from the samples to obtain posterior distributions for the ratios.
#     The ratios and their intervals are then visualized across disease categories.

# Prepare the data variables for analysis.

IG_all <- antibody_family_complete_Ig[,1:4]
IG_all$Total <- antibody_family_complete_Ig$Total
IG_all$Family <- as.matrix(antibody_family_complete_Ig[,5:16])
IG_all$Ig <- factor(IG_all$Ig, levels = c('Presort','IgAP','IgGP'))

# Conduct the hierarchical multinomial model, and print the results summary.

family_all_mod <- brm(Family|trials(Total) ~ 1 + Classification_CE*Stage*Ig+(1|Subject), family = multinomial(), data = IG_all_deID, cores = 4, prior = set_prior('normal(0,5)', class = 'b'), control = list(max_treedepth = 15))
summary(family_IGA_mod)

write.csv(IG_all_deID, file = 'antibody_family_data_AW_deID.csv', row.names = FALSE)

# Prepare the posterior predictions and generate an interval plot of the population proportional abundances for each taxon.

new_fitted_data <- data.frame(Classification_CE = character(8), Stage = character(8), Ig = character(8), Total = numeric(8))
new_fitted_data$Classification_CE <- as.factor(rep(c('DRE','ARE','IRE','Healthy'),2))
new_fitted_data$Stage <- as.factor(rep(c(rep('Before',4),rep('After',4)),1))
new_fitted_data$Ig <- as.factor(rep('Presort',8))
new_fitted_data$Total <- 1
new_fitted_obs1 <- fitted(family_all_mod, newdata = new_fitted_data, re_formula = NA, probs = c(0.05,0.95))

new_fitted_panelnames <- dimnames(new_fitted_obs1)[3]
new_fitted_colnames <- colnames(new_fitted_obs1)

new_pred_presort <- as.data.frame(rbind(new_fitted_obs1[,,1],new_fitted_obs1[,,2],new_fitted_obs1[,,3],new_fitted_obs1[,,4],new_fitted_obs1[,,5],new_fitted_obs1[,,6],new_fitted_obs1[,,7],new_fitted_obs1[,,8],new_fitted_obs1[,,9],new_fitted_obs1[,,10],new_fitted_obs1[,,11],new_fitted_obs1[,,12]))
colnames(new_pred_presort) <- new_fitted_colnames

new_pred_presort$Classification_CE <- rep(new_fitted_data$Classification_CE,12)
new_pred_presort$Stage <- rep(new_fitted_data$Stage,12)
new_pred_presort$Ig <- rep(new_fitted_data$Ig,12)

taxa_names <- c(rep(substr(new_fitted_panelnames[[1]][1],7,nchar(new_fitted_panelnames[[1]][1])-1),8),
                rep(substr(new_fitted_panelnames[[1]][2],7,nchar(new_fitted_panelnames[[1]][2])-1),8),
                rep(substr(new_fitted_panelnames[[1]][3],7,nchar(new_fitted_panelnames[[1]][3])-1),8),
                rep(substr(new_fitted_panelnames[[1]][4],7,nchar(new_fitted_panelnames[[1]][4])-1),8),
                rep(substr(new_fitted_panelnames[[1]][5],7,nchar(new_fitted_panelnames[[1]][5])-1),8),
                rep(substr(new_fitted_panelnames[[1]][6],7,nchar(new_fitted_panelnames[[1]][6])-1),8),
                rep(substr(new_fitted_panelnames[[1]][7],7,nchar(new_fitted_panelnames[[1]][7])-1),8),
                rep(substr(new_fitted_panelnames[[1]][8],7,nchar(new_fitted_panelnames[[1]][8])-1),8),
                rep(substr(new_fitted_panelnames[[1]][9],7,nchar(new_fitted_panelnames[[1]][9])-1),8),
                rep(substr(new_fitted_panelnames[[1]][10],7,nchar(new_fitted_panelnames[[1]][10])-1),8),
                rep(substr(new_fitted_panelnames[[1]][11],7,nchar(new_fitted_panelnames[[1]][11])-1),8),
                rep(substr(new_fitted_panelnames[[1]][12],7,nchar(new_fitted_panelnames[[1]][12])-1),8))

new_pred_presort$Taxon <- taxa_names

new_pred_presort$Stage <- factor(new_pred_presort$Stage, levels = c('Before','After'))
new_pred_presort$Classification_CE <- factor(new_pred_presort$Classification_CE, levels = c('Healthy','ARE','DRE','IRE'))

taxon_abundance_plot <- ggplot(data = new_pred_presort, aes(x = Classification_CE, y = Estimate, color = Stage)) + geom_point(position = position_dodge(width = 0.5)) + geom_errorbar(aes(ymin = Q5, ymax = Q95), position = position_dodge(width = 0.5)) + facet_wrap(~Taxon) + ylab('Posterior Proportional Abundance') + theme(axis.title.x=element_blank()) + theme(strip.text = element_text(face = "italic"))
ggsave(filename = 'taxon_abundance_plot.tif', plot = taxon_abundance_plot, device = 'tiff', dpi = 100, width = 7, height = 7)

print(new_pred_presort)

# Construct a matrix of the samples for the population proportional abundances for both the IgA and presort conditions, and then plot these posterior distributions.

new_fitted_IgA_data <- new_fitted_data
new_fitted_IgA_data$Ig <- 'IgAP'

IgAP_samples <- fitted(family_all_mod, newdata = new_fitted_IgA_data, summary = FALSE, re_formula = NA)
IgAN_samples <- fitted(family_all_mod, newdata = new_fitted_data, summary = FALSE, re_formula = NA)
IgA_ratio <- IgAP_samples/IgAN_samples

IgA_ratio_vec1 <- as.data.frame(c(IgA_ratio[,1,1],IgA_ratio[,2,1],IgA_ratio[,3,1],IgA_ratio[,4,1],IgA_ratio[,5,1],IgA_ratio[,6,1],IgA_ratio[,7,1],IgA_ratio[,8,1]))
colnames(IgA_ratio_vec1) <- 'sampled_value'
IgA_ratio_vec1$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec1$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec1$Taxon <- rep(substr(new_fitted_panelnames[[1]][1],7,nchar(new_fitted_panelnames[[1]][1])-1),32000)

IgA_ratio_vec2 <- as.data.frame(c(IgA_ratio[,1,2],IgA_ratio[,2,2],IgA_ratio[,3,2],IgA_ratio[,4,2],IgA_ratio[,5,2],IgA_ratio[,6,2],IgA_ratio[,7,2],IgA_ratio[,8,2]))
colnames(IgA_ratio_vec2) <- 'sampled_value'
IgA_ratio_vec2$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec2$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec2$Taxon <- rep(substr(new_fitted_panelnames[[1]][2],7,nchar(new_fitted_panelnames[[1]][2])-1),32000)

IgA_ratio_vec3 <- as.data.frame(c(IgA_ratio[,1,3],IgA_ratio[,2,3],IgA_ratio[,3,3],IgA_ratio[,4,3],IgA_ratio[,5,3],IgA_ratio[,6,3],IgA_ratio[,7,3],IgA_ratio[,8,3]))
colnames(IgA_ratio_vec3) <- 'sampled_value'
IgA_ratio_vec3$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec3$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec3$Taxon <- rep(substr(new_fitted_panelnames[[1]][3],7,nchar(new_fitted_panelnames[[1]][3])-1),32000)

IgA_ratio_vec4 <- as.data.frame(c(IgA_ratio[,1,4],IgA_ratio[,2,4],IgA_ratio[,3,4],IgA_ratio[,4,4],IgA_ratio[,5,4],IgA_ratio[,6,4],IgA_ratio[,7,4],IgA_ratio[,8,4]))
colnames(IgA_ratio_vec4) <- 'sampled_value'
IgA_ratio_vec4$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec4$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec4$Taxon <- rep(substr(new_fitted_panelnames[[1]][4],7,nchar(new_fitted_panelnames[[1]][4])-1),32000)

IgA_ratio_vec5 <- as.data.frame(c(IgA_ratio[,1,5],IgA_ratio[,2,5],IgA_ratio[,3,5],IgA_ratio[,4,5],IgA_ratio[,5,5],IgA_ratio[,6,5],IgA_ratio[,7,5],IgA_ratio[,8,5]))
colnames(IgA_ratio_vec5) <- 'sampled_value'
IgA_ratio_vec5$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec5$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec5$Taxon <- rep(substr(new_fitted_panelnames[[1]][5],7,nchar(new_fitted_panelnames[[1]][5])-1),32000)

IgA_ratio_vec6 <- as.data.frame(c(IgA_ratio[,1,6],IgA_ratio[,2,6],IgA_ratio[,3,6],IgA_ratio[,4,6],IgA_ratio[,5,6],IgA_ratio[,6,6],IgA_ratio[,7,6],IgA_ratio[,8,6]))
colnames(IgA_ratio_vec6) <- 'sampled_value'
IgA_ratio_vec6$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec6$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec6$Taxon <- rep(substr(new_fitted_panelnames[[1]][6],7,nchar(new_fitted_panelnames[[1]][6])-1),32000)

IgA_ratio_vec7 <- as.data.frame(c(IgA_ratio[,1,7],IgA_ratio[,2,7],IgA_ratio[,3,7],IgA_ratio[,4,7],IgA_ratio[,5,7],IgA_ratio[,6,7],IgA_ratio[,7,7],IgA_ratio[,8,7]))
colnames(IgA_ratio_vec7) <- 'sampled_value'
IgA_ratio_vec7$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec7$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec7$Taxon <- rep(substr(new_fitted_panelnames[[1]][7],7,nchar(new_fitted_panelnames[[1]][7])-1),32000)

IgA_ratio_vec8 <- as.data.frame(c(IgA_ratio[,1,8],IgA_ratio[,2,8],IgA_ratio[,3,8],IgA_ratio[,4,8],IgA_ratio[,5,8],IgA_ratio[,6,8],IgA_ratio[,7,8],IgA_ratio[,8,8]))
colnames(IgA_ratio_vec8) <- 'sampled_value'
IgA_ratio_vec8$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec8$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec8$Taxon <- rep(substr(new_fitted_panelnames[[1]][8],7,nchar(new_fitted_panelnames[[1]][8])-1),32000)

IgA_ratio_vec9 <- as.data.frame(c(IgA_ratio[,1,9],IgA_ratio[,2,9],IgA_ratio[,3,9],IgA_ratio[,4,9],IgA_ratio[,5,9],IgA_ratio[,6,9],IgA_ratio[,7,9],IgA_ratio[,8,9]))
colnames(IgA_ratio_vec9) <- 'sampled_value'
IgA_ratio_vec9$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec9$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec9$Taxon <- rep(substr(new_fitted_panelnames[[1]][9],7,nchar(new_fitted_panelnames[[1]][9])-1),32000)

IgA_ratio_vec10 <- as.data.frame(c(IgA_ratio[,1,10],IgA_ratio[,2,10],IgA_ratio[,3,10],IgA_ratio[,4,10],IgA_ratio[,5,10],IgA_ratio[,6,10],IgA_ratio[,7,10],IgA_ratio[,8,10]))
colnames(IgA_ratio_vec10) <- 'sampled_value'
IgA_ratio_vec10$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec10$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec10$Taxon <- rep(substr(new_fitted_panelnames[[1]][10],7,nchar(new_fitted_panelnames[[1]][10])-1),32000)

IgA_ratio_vec11 <- as.data.frame(c(IgA_ratio[,1,11],IgA_ratio[,2,11],IgA_ratio[,3,11],IgA_ratio[,4,11],IgA_ratio[,5,11],IgA_ratio[,6,11],IgA_ratio[,7,11],IgA_ratio[,8,11]))
colnames(IgA_ratio_vec11) <- 'sampled_value'
IgA_ratio_vec11$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec11$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec11$Taxon <- rep(substr(new_fitted_panelnames[[1]][11],7,nchar(new_fitted_panelnames[[1]][11])-1),32000)

IgA_ratio_vec12 <- as.data.frame(c(IgA_ratio[,1,12],IgA_ratio[,2,12],IgA_ratio[,3,12],IgA_ratio[,4,12],IgA_ratio[,5,12],IgA_ratio[,6,12],IgA_ratio[,7,12],IgA_ratio[,8,12]))
colnames(IgA_ratio_vec12) <- 'sampled_value'
IgA_ratio_vec12$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgA_ratio_vec12$Stage <- c(rep('Before',16000),rep('After',16000))
IgA_ratio_vec12$Taxon <- rep(substr(new_fitted_panelnames[[1]][12],7,nchar(new_fitted_panelnames[[1]][12])-1),32000)

IgA_total_vec <- rbind(IgA_ratio_vec1,IgA_ratio_vec2,IgA_ratio_vec3,IgA_ratio_vec4,IgA_ratio_vec5,IgA_ratio_vec6,IgA_ratio_vec7,IgA_ratio_vec8,IgA_ratio_vec9,IgA_ratio_vec10,IgA_ratio_vec11,IgA_ratio_vec12)

IgA_total_estimates <- IgA_total_vec

IgA_total_estimates$Classification_CE <- factor(IgA_total_estimates$Classification_CE, levels = c('IRE', 'DRE', 'ARE', 'Healthy'), ordered = TRUE)
IgA_total_estimates$Stage <- factor(IgA_total_estimates$Stage, levels = c('Before','After'))

IgA_ratio_plot <- ggplot(data = IgA_total_estimates, aes(x = sampled_value, y = Classification_CE, color = Stage)) + stat_interval(.width = 0.5, position = position_dodge2(width = 0.6, reverse = TRUE)) + stat_interval(.width = 0.95, position = position_dodge2(width = 0.6, reverse = TRUE), alpha = 0.4) +  scale_x_log10(breaks = c(0.01,1,100), labels = c(0.01,1,100)) + coord_cartesian(xlim = c(0.001,1000)) + facet_wrap(~Taxon) + geom_vline(xintercept = 1) + xlab('Immunoglobulin-A Enrichment Ratio') + ylab('') + theme(strip.text = element_text(face = "italic"), legend.position = 'none')

# Construct a matrix of the samples for the population proportional abundances for both the IgG and presort conditions, and then plot these posterior distributions.

new_fitted_IgG_data <- new_fitted_data
new_fitted_IgG_data$Ig <- 'IgGP'

IgGP_samples <- fitted(family_all_mod, newdata = new_fitted_IgG_data, summary = FALSE, re_formula = NA)
IgGN_samples <- fitted(family_all_mod, newdata = new_fitted_data, summary = FALSE, re_formula = NA)
IgG_ratio <- IgGP_samples/IgGN_samples


IgG_ratio_vec1 <- as.data.frame(c(IgG_ratio[,1,1],IgG_ratio[,2,1],IgG_ratio[,3,1],IgG_ratio[,4,1],IgG_ratio[,5,1],IgG_ratio[,6,1],IgG_ratio[,7,1],IgG_ratio[,8,1]))
colnames(IgG_ratio_vec1) <- 'sampled_value'
IgG_ratio_vec1$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec1$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec1$Taxon <- rep(substr(new_fitted_panelnames[[1]][1],7,nchar(new_fitted_panelnames[[1]][1])-1),32000)

IgG_ratio_vec2 <- as.data.frame(c(IgG_ratio[,1,2],IgG_ratio[,2,2],IgG_ratio[,3,2],IgG_ratio[,4,2],IgG_ratio[,5,2],IgG_ratio[,6,2],IgG_ratio[,7,2],IgG_ratio[,8,2]))
colnames(IgG_ratio_vec2) <- 'sampled_value'
IgG_ratio_vec2$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec2$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec2$Taxon <- rep(substr(new_fitted_panelnames[[1]][2],7,nchar(new_fitted_panelnames[[1]][2])-1),32000)

IgG_ratio_vec3 <- as.data.frame(c(IgG_ratio[,1,3],IgG_ratio[,2,3],IgG_ratio[,3,3],IgG_ratio[,4,3],IgG_ratio[,5,3],IgG_ratio[,6,3],IgG_ratio[,7,3],IgG_ratio[,8,3]))
colnames(IgG_ratio_vec3) <- 'sampled_value'
IgG_ratio_vec3$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec3$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec3$Taxon <- rep(substr(new_fitted_panelnames[[1]][3],7,nchar(new_fitted_panelnames[[1]][3])-1),32000)

IgG_ratio_vec4 <- as.data.frame(c(IgG_ratio[,1,4],IgG_ratio[,2,4],IgG_ratio[,3,4],IgG_ratio[,4,4],IgG_ratio[,5,4],IgG_ratio[,6,4],IgG_ratio[,7,4],IgG_ratio[,8,4]))
colnames(IgG_ratio_vec4) <- 'sampled_value'
IgG_ratio_vec4$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec4$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec4$Taxon <- rep(substr(new_fitted_panelnames[[1]][4],7,nchar(new_fitted_panelnames[[1]][4])-1),32000)

IgG_ratio_vec5 <- as.data.frame(c(IgG_ratio[,1,5],IgG_ratio[,2,5],IgG_ratio[,3,5],IgG_ratio[,4,5],IgG_ratio[,5,5],IgG_ratio[,6,5],IgG_ratio[,7,5],IgG_ratio[,8,5]))
colnames(IgG_ratio_vec5) <- 'sampled_value'
IgG_ratio_vec5$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec5$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec5$Taxon <- rep(substr(new_fitted_panelnames[[1]][5],7,nchar(new_fitted_panelnames[[1]][5])-1),32000)

IgG_ratio_vec6 <- as.data.frame(c(IgG_ratio[,1,6],IgG_ratio[,2,6],IgG_ratio[,3,6],IgG_ratio[,4,6],IgG_ratio[,5,6],IgG_ratio[,6,6],IgG_ratio[,7,6],IgG_ratio[,8,6]))
colnames(IgG_ratio_vec6) <- 'sampled_value'
IgG_ratio_vec6$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec6$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec6$Taxon <- rep(substr(new_fitted_panelnames[[1]][6],7,nchar(new_fitted_panelnames[[1]][6])-1),32000)

IgG_ratio_vec7 <- as.data.frame(c(IgG_ratio[,1,7],IgG_ratio[,2,7],IgG_ratio[,3,7],IgG_ratio[,4,7],IgG_ratio[,5,7],IgG_ratio[,6,7],IgG_ratio[,7,7],IgG_ratio[,8,7]))
colnames(IgG_ratio_vec7) <- 'sampled_value'
IgG_ratio_vec7$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec7$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec7$Taxon <- rep(substr(new_fitted_panelnames[[1]][7],7,nchar(new_fitted_panelnames[[1]][7])-1),32000)

IgG_ratio_vec8 <- as.data.frame(c(IgG_ratio[,1,8],IgG_ratio[,2,8],IgG_ratio[,3,8],IgG_ratio[,4,8],IgG_ratio[,5,8],IgG_ratio[,6,8],IgG_ratio[,7,8],IgG_ratio[,8,8]))
colnames(IgG_ratio_vec8) <- 'sampled_value'
IgG_ratio_vec8$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec8$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec8$Taxon <- rep(substr(new_fitted_panelnames[[1]][8],7,nchar(new_fitted_panelnames[[1]][8])-1),32000)

IgG_ratio_vec9 <- as.data.frame(c(IgG_ratio[,1,9],IgG_ratio[,2,9],IgG_ratio[,3,9],IgG_ratio[,4,9],IgG_ratio[,5,9],IgG_ratio[,6,9],IgG_ratio[,7,9],IgG_ratio[,8,9]))
colnames(IgG_ratio_vec9) <- 'sampled_value'
IgG_ratio_vec9$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec9$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec9$Taxon <- rep(substr(new_fitted_panelnames[[1]][9],7,nchar(new_fitted_panelnames[[1]][9])-1),32000)

IgG_ratio_vec10 <- as.data.frame(c(IgG_ratio[,1,10],IgG_ratio[,2,10],IgG_ratio[,3,10],IgG_ratio[,4,10],IgG_ratio[,5,10],IgG_ratio[,6,10],IgG_ratio[,7,10],IgG_ratio[,8,10]))
colnames(IgG_ratio_vec10) <- 'sampled_value'
IgG_ratio_vec10$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec10$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec10$Taxon <- rep(substr(new_fitted_panelnames[[1]][10],7,nchar(new_fitted_panelnames[[1]][10])-1),32000)

IgG_ratio_vec11 <- as.data.frame(c(IgG_ratio[,1,11],IgG_ratio[,2,11],IgG_ratio[,3,11],IgG_ratio[,4,11],IgG_ratio[,5,11],IgG_ratio[,6,11],IgG_ratio[,7,11],IgG_ratio[,8,11]))
colnames(IgG_ratio_vec11) <- 'sampled_value'
IgG_ratio_vec11$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec11$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec11$Taxon <- rep(substr(new_fitted_panelnames[[1]][11],7,nchar(new_fitted_panelnames[[1]][11])-1),32000)

IgG_ratio_vec12 <- as.data.frame(c(IgG_ratio[,1,12],IgG_ratio[,2,12],IgG_ratio[,3,12],IgG_ratio[,4,12],IgG_ratio[,5,12],IgG_ratio[,6,12],IgG_ratio[,7,12],IgG_ratio[,8,12]))
colnames(IgG_ratio_vec12) <- 'sampled_value'
IgG_ratio_vec12$Classification_CE <- c(rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000),rep('DRE',4000),rep('ARE',4000),rep('IRE',4000),rep('Healthy',4000))
IgG_ratio_vec12$Stage <- c(rep('Before',16000),rep('After',16000))
IgG_ratio_vec12$Taxon <- rep(substr(new_fitted_panelnames[[1]][12],7,nchar(new_fitted_panelnames[[1]][12])-1),32000)

IgG_total_vec <- rbind(IgG_ratio_vec1,IgG_ratio_vec2,IgG_ratio_vec3,IgG_ratio_vec4,IgG_ratio_vec5,IgG_ratio_vec6,IgG_ratio_vec7,IgG_ratio_vec8,IgG_ratio_vec9,IgG_ratio_vec10,IgG_ratio_vec11,IgG_ratio_vec12)

IgG_total_estimates <- IgG_total_vec

IgG_total_estimates$Classification_CE <- factor(IgG_total_estimates$Classification_CE, levels = c('IRE', 'DRE', 'ARE', 'Healthy'), ordered = TRUE)
IgG_total_estimates$Stage <- factor(IgG_total_estimates$Stage, levels = c('Before','After'))

IgG_ratio_plot <- ggplot(data = IgG_total_estimates, aes(x = sampled_value, y = Classification_CE, color = Stage)) + stat_interval(.width = 0.5, position = position_dodge2(width = 0.6, reverse = TRUE)) + stat_interval(.width = 0.95, position = position_dodge2(width = 0.6, reverse = TRUE), alpha = 0.4) +  scale_x_log10(breaks = c(0.01,1,100), labels = c(0.01,1,100)) + coord_cartesian(xlim = c(0.001,1000)) + facet_wrap(~Taxon) + geom_vline(xintercept = 1) + xlab('Immunoglobulin-G Enrichment Ratio') + ylab('') + theme(strip.text = element_text(face = "italic"), legend.position = 'none')

#

Ig_joined_plots <- ggarrange(IgA_ratio_plot, IgG_ratio_plot, ncol = 1, common.legend = TRUE, legend = 'right')
ggsave(filename = 'Ig_joined_plots.tif', plot = Ig_joined_plots, device = 'tiff', dpi = 100, width = 7, height = 7)

# Prepare the posterior predictions for each subject, only for the conditions actually observed in the experiment. This is equivalent to the dataset (but without duplicate rows). 

subject_predictors_all <- unique(antibody_family_complete_Ig[antibody_family_complete_Ig$Ig == 'Presort',1:4])
subject_predictors_all$Classification_CE <- factor(subject_predictors_all$Classification_CE, levels = c('Healthy','ARE','DRE','IRE'), ordered = FALSE)
subject_predictors_all <- subject_predictors_all[order(subject_predictors_all$Classification_CE, subject_predictors_all$Dog),]
subject_predictors_all$Total <- 1

subject_all_pred <- fitted(family_all_mod, newdata = subject_predictors_all)

new_subject_presort <- as.data.frame(rbind(subject_all_pred[,,1],subject_all_pred[,,2],subject_all_pred[,,3],subject_all_pred[,,4],subject_all_pred[,,5],subject_all_pred[,,6],subject_all_pred[,,7],subject_all_pred[,,8],subject_all_pred[,,9],subject_all_pred[,,10],subject_all_pred[,,11],subject_all_pred[,,12]))
colnames(new_subject_presort) <- new_fitted_colnames

new_subject_presort$Stage <- factor(rep(subject_predictors_all$Stage,12), levels = c('Before','After'))
new_subject_presort$Ig <- factor(rep(subject_predictors_all$Ig,12))
new_subject_presort$Dog <- factor(rep(subject_predictors_all$Dog,12), levels = unique(subject_predictors_all$Dog))
new_subject_presort$Classification_CE <- factor(rep(subject_predictors_all$Classification_CE,12), levels = c('Healthy','ARE','DRE','IRE'))

subject_taxa_names <- c(rep(substr(new_fitted_panelnames[[1]][1],7,nchar(new_fitted_panelnames[[1]][1])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][2],7,nchar(new_fitted_panelnames[[1]][2])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][3],7,nchar(new_fitted_panelnames[[1]][3])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][4],7,nchar(new_fitted_panelnames[[1]][4])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][5],7,nchar(new_fitted_panelnames[[1]][5])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][6],7,nchar(new_fitted_panelnames[[1]][6])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][7],7,nchar(new_fitted_panelnames[[1]][7])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][8],7,nchar(new_fitted_panelnames[[1]][8])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][9],7,nchar(new_fitted_panelnames[[1]][9])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][10],7,nchar(new_fitted_panelnames[[1]][10])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][11],7,nchar(new_fitted_panelnames[[1]][11])-1),54),
                        rep(substr(new_fitted_panelnames[[1]][12],7,nchar(new_fitted_panelnames[[1]][12])-1),54))

new_subject_presort$Taxon <- subject_taxa_names

subject_estimates_plot <- ggplot(data = new_subject_presort, aes(x = Dog, y = Estimate, ymin = Q5, ymax = Q95, color = Classification_CE, shape = Stage)) +  geom_pointrange() + facet_wrap(~Taxon) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab('Posterior Proportional Abundance') + theme(strip.text = element_text(face = "italic"))
ggsave(filename = 'subject_estimates_plot.tif', plot = subject_estimates_plot, device = 'tiff', dpi = 100, width = 8, height = 8)
