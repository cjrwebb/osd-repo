###############################################################################
#########      ANALYSIS SCRIPT FOR: FUNDING AND DEPRIVATION AS     ############
#########  DETERMINANTS OF CHILD WELFARE AND PROTECTION SERVICES   ############
######### QUALITY: A MULTILEVEL ANALYSIS OF EVALUATIONS IN ENGLAND ############
############################################################################### 
#########     CALUM J R WEBB, DAVARA L BENNETT, PAUL W B BYWATERS   ###########
###############################################################################

# load important packages for data manipulation, model estimation, and marginal
# effects estimation.
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(lavaan)
library(psych)
library(performance)
library(caret)

# read in data of 374 Ofsted quality inspections linked to local authority 
# spending and deprivation data
analysis_data <- read_csv(file = "data.csv")

# Create null model (random effects only)
mod_null <- glmmTMB(data = analysis_data, 
                    formula = good_or_outstanding ~ (1 | UTLA18CD) + (1 | year),
                    family = binomial(link = "logit"))

# Create model with untransformed predictors (required for ggeffects marginal effects)
mod1 <- glmmTMB(data =    analysis_data, 
                formula = good_or_outstanding ~ exp_ehfs_perchild + exp_sg_perchild + 
                          imd_decile + (1 | UTLA18CD) + (1 | year),
                family =  binomial(link = "logit")) 

# Create model with centred expenditure in £100s
mod2 <- glmmTMB(data = analysis_data, 
                formula = good_or_outstanding ~ I(scale(exp_ehfs_perchild/100, scale = FALSE)) + 
                          I(scale(exp_sg_perchild/100, scale = FALSE)) + imd_decile + 
                          (1 | UTLA18CD) + (1 | year),
                family = binomial(link = "logit")) 

# Print model output
summary(mod_null)
summary(mod1)
summary(mod2)

# Likelihood ratio test between null model and model using spend and deprivation
performance_lrt(mod_null, mod1)

################### Calculate accuracy of model and confusion matrix
#' Calculate the predicted probabilities for each inspection based on
#' model parameters.
predicted_percent <- exp(predict(mod1))/(1+exp(predict(mod1)))

#' Combine actual outcomes with predicted outcomes using a cutoff of 0.5
#' (predictions greater than or equal to 50% = Good/Outstanding, predictions
#' less than 50% = Requires Improvement/Inadequate).
conf_matrix_data <- tibble(actual = analysis_data$good_or_outstanding,
                           model = predicted_percent) %>%
  mutate(model = ifelse(model < 0.5, 0, 1)) %>%
  mutate_all(~as.factor(ifelse(. == 1, "Good/Outstanding", "Inadequate/RI")))

confusionMatrix(conf_matrix_data$model, reference = conf_matrix_data$actual)

# Plot predicted probabilities and save as .eps for each predictor
(plot1 <- ggpredict(mod1, terms = c("exp_ehfs_perchild[100:1000]")) %>% 
    plot(limits = c(0,1), ci.style = "dash") + 
    ylab("Probability of Good or Outstanding Ofsted Inspection") +
    xlab("Expenditure on Preventative Services per Child (£)\n(Range: Min & Max in Sample, rounded)") +
    ggtitle("")
)

ggsave(plot = plot1, filename = "plot_1_ehfs.eps", width = 6, height = 4)

(plot2 <- ggpredict(mod1, terms = c("exp_sg_perchild[30:650]")) %>% 
    plot(limits = c(0,1), ci.style = "dash") + 
    ylab("Probability of Good or Outstanding Ofsted Inspection") +
    xlab("Expenditure on Safeguarding per Child (£)\n(Range: Min & Max in Sample, rounded)") +
    ggtitle(""))

ggsave(plot = plot2, filename = "plot_2_sg.eps", width = 6, height = 4)

(plot3 <- ggpredict(mod1, terms = c("imd_decile[1:10]")) %>% 
    plot(limits = c(0,1), ci.style = "dash") + 
    scale_x_continuous(breaks = seq(1, 10, 1)) +
    ylab("Probability of Good or Outstanding Ofsted Inspection") +
    xlab("Decile of Deprivation (Measured by IMD 2019)\nHigher = More Deprived") +
    ggtitle(""))

ggsave(plot = plot3, filename = "plot_3_imd.eps", width = 6, height = 4)

############### Internal reliability and communality across frameworks and areas of 
############### evaluation

#' Read in and recode ofsted subdomain/areas of inspection data for all
#' frameworks. Drop NA cases.

# Including ILACS (limits N = 70)
corres_data <- read_csv("ofsted_subdomains.csv") %>%
  mutate_at(vars(slac_sg:sif_mng), ~case_when(. == "Outstanding" ~ 3,
                                              . == "Good" ~ 2,
                                              . == "Requires improvement to be good" ~ 1,
                                              . == "Inadequate" ~ 0)) %>%
  select(-ilacs_oe:-ilacs_cla) %>%
  drop_na() %>%
  select(-id)

# Excluding ILACS (N = 150)
corres_data_ilacs <- read_csv("ofsted_subdomains.csv") %>%
  mutate_at(vars(slac_sg:ilacs_cla), ~case_when(. == "Outstanding" ~ 3,
                                                . == "Good" ~ 2,
                                                . == "Requires improvement to be good" ~ 1,
                                                . == "Inadequate" ~ 0)) %>%
  drop_na() %>%
  select(-id)

# Alpha and Gutmann's Lambda-6 for all frameworks/areas of evaluation
alpha(corres_data_ilacs)

# Alpha and Gutmann's Lambda-6 for SLAC - SIF
alpha(corres_data)

# Alpha and Gutmann's Lambda-6 for SIF
psych::alpha(corres_data[, 3:8])

# Alpha and Gutmann's Lambda-6 for SLAC
psych::alpha(corres_data[, 1:2])

# Alpha and Gutmann's Lambda-6 for ILACS
psych::alpha(corres_data_ilacs[, 9:12])


######################### Measures of Communality using 1F CFA

# SLAC Spearman's rho and Wilcox test for significance
cor(corres_data$slac_sg, corres_data$slac_cla, method = "spearman")
wilcox.test(corres_data$slac_sg, corres_data$slac_cla)

# Define CFA for SIF domains
cfa_1_model <- "
sif =~ 1*sif_oe + sif_cin + sif_cla + sif_adp + sif_cl + 1*sif_mng
"

# Define CFA for ILACS domains
cfa_2_model <- "
ilacs =~ ilacs_oe + ilacs_ldr + ilacs_cin + ilacs_cla
ilacs_oe ~~ 0*ilacs_oe
"

# run and print output of one factor CFA models for SIF and ILACS
sem(model = cfa_1_model, data = corres_data, 
    ordered = c("sif_oe", "sif_cin", "sif_cla", "sif_adp", "sif_cl", "sif_mng")) %>% 
  summary(., fit.measures = TRUE, standardized = TRUE)

sem(model = cfa_2_model, data = corres_data_ilacs, 
    ordered = c("ilacs_oe", "ilacs_ldr", "ilacs_cin", "ilacs_cla")) %>% 
  summary(., fit.measures = TRUE, standardized = TRUE)
