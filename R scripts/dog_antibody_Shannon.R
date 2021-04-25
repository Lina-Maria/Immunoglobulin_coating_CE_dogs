# AP. Woodward and L. Martinez-Lopez, 2020-21, University of Melbourne.
# Conducts analysis of alpha-diversity (Shannon index) from the microbiome dataset including antibody status, disease classification, and stage
#     as predictors.
# Then plots and summarizes the results.

## Load packages
library(lme4)
library(emmeans)
library(MuMIn)

# Conduct a linear mixed model for the log(Shannon) with random intercept for Dog.
#     log-transformation is natural for Shannon index as it is positive-nonzero.
#     The transformation is readily justified by the diagnostic plots.
# The full set of possible interacting terms are used here.

lmer(log(Shannon)~1+Ig*Classification_CE*Stage+(1|Dog), data = ShannonIg)

Shannon_mod <-lmer(log(Shannon)~1+Ig*Classification_CE*Stage+(1|Dog), data = ShannonIg)

summary(Shannon_mod)
confint(Shannon_mod)

hist(residuals(Shannon_mod))
plot(fitted(Shannon_mod), residuals(Shannon_mod))

r.squaredGLMM(Shannon_mod)

# Using the 'emmeans' package, develop predicted expected value of Shannon under each condition, and then generate an interval plot of the results.
#     Also print the 'emmeans' table.

Shannon_post_means <- emmeans(Shannon_mod, specs = c('Ig','Classification_CE','Stage'), type = 'response')

Shannon_post_data <- as.data.frame(Shannon_post_means)
Shannon_post_data$Classification_CE <- factor(Shannon_post_data$Classification_CE, levels = c('Healthy','ARE','FRE','IRE'))
Shannon_post_data$Stage <- factor(Shannon_post_data$Stage, levels = c('Before','After'))
Shannon_post_data$Ig <- factor(Shannon_post_data$Ig, levels = c('Input','IgAp','IgGp'))

ggplot(data = Shannon_post_data, aes(x = Classification_CE, color = Stage, y = response, ymin = lower.CL, ymax = upper.CL)) + geom_point(shape = 2, position = position_dodge(width = 0.7)) + geom_errorbar(position = position_dodge(width = 0.7)) + facet_wrap(~Ig, nrow = 1) + theme(axis.title.x=element_blank()) + ylab('Estimated Shannon Index') + ylim(c(1.5,4.5))
print(Shannon_post_means)

#Ig    Classification_CE Stage  response    SE    df lower.CL upper.CL
#IgAp  ARE               After      2.05 0.192 106.2     1.70     2.47
#IgGp  ARE               After      2.03 0.191 106.2     1.68     2.45
#Input ARE               After      2.86 0.226  59.2     2.45     3.35
#IgAp  FRE               After      2.45 0.216 104.7     2.06     2.92
#IgGp  FRE               After      2.41 0.212 104.7     2.03     2.87
#Input FRE               After      3.17 0.234  58.3     2.73     3.67
#IgAp  Healthy           After      2.29 0.172  87.0     1.98     2.66
#IgGp  Healthy           After      2.05 0.161 102.5     1.75     2.39
#Input Healthy           After      2.82 0.181  50.8     2.48     3.21
#IgAp  IRE               After      2.49 0.422 134.0     1.78     3.48
#IgGp  IRE               After      2.99 0.508 134.0     2.14     4.19
#Input IRE               After      3.31 0.464  75.7     2.50     4.38
#IgAp  ARE               Before     2.66 0.229  80.4     2.24     3.16
#IgGp  ARE               Before     2.57 0.221  80.4     2.16     3.05
#Input ARE               Before     3.07 0.229  48.5     2.64     3.56
#IgAp  FRE               Before     2.32 0.176  65.2     1.99     2.70
#IgGp  FRE               Before     2.41 0.183  65.2     2.07     2.80
#Input FRE               Before     2.70 0.180  40.5     2.35     3.09
#IgAp  Healthy           Before     1.87 0.148 102.5     1.60     2.19
#IgGp  Healthy           Before     2.13 0.162  93.5     1.83     2.47
#Input Healthy           Before     2.58 0.168  54.4     2.26     2.94
#IgAp  IRE               Before     2.09 0.296  78.1     1.58     2.77
#IgGp  IRE               Before     2.16 0.306  78.1     1.63     2.87
#Input IRE               Before     2.78 0.345  49.0     2.16     3.56

