library(tidyverse)
library(tableone)
library(patchwork)
library(brms)
library(tidybayes)
library(bayesplot)
library(future)
library(corrplot)
library(lme4)

##### Load data #####
df1 <- readRDS("Data/REVEALS_formatted_data.rds")


##### Data preparation #####
df1$uin <- as.factor(df1$uin)
# extract unique IDS
id_list <- unique(df1$uin)

# Make a baseline dataframe with first visit of each individual
df_base <- df1[ df1$ax_number == 1, ]

# Determine time of last outcome per individual
last_outcome <- df1 %>% 
    group_by(uin) %>%
    filter(ax_date == max(ax_date)) %>%
    ungroup()
df_base$follow_up_time <- last_outcome[ match(df_base$uin, last_outcome$uin), ]$days_from_baseline

# Determine number of longitudinal measurements per person
outcome_counts <- df1 %>% count(uin)
df_base$count_long <- outcome_counts[ match(df_base$uin, outcome_counts$uin), ]$n

# Exclude any individuals not at stage 2 or 3 at presentation
df_base <- df_base[ df_base$staging %in% c("Kings 2", "Kings 3") | is.na(df_base$staging), ]
df1 <- df1[df1$uin %in% df_base$uin, ]

# Re-factor the staging variable
df_base$staging <- as.factor( as.character(df_base$staging))

# Rename outcome variables with shorter names to make model specification more concise
df1 <- df1 %>% rename(out_fvc = fvc_max_l,
                      out_svc = svc_l,
                      out_snip = snip_max_score_cm_h2o,
                      out_peak = peak_cough_flow_max_score)


######## Comparisons of respiratory outcomes using REVEALS dataset #######

#### Make Table 1 of descriptive statistics ####
vars <- c("gender", "diagnosis", "staging","height_cm", 
          "weight_kg", "age_dx", "follow_up_time", "count_long")
non_par <- c("follow_up_time", "count_long")

tab1 <- CreateTableOne(vars, strata = "cohort", data = df_base )
tab1b <- print(tab1, nonnormal = non_par, exact = c("gender", "diagnosis", "staging"),
               quote = FALSE, noSpaces = TRUE,
               printToggle = FALSE, showAllLevels = TRUE)
tab1_total <- CreateTableOne(vars, data = df_base)
tab1_total <- print(tab1_total, nonnormal = non_par, exact = c("gender", "diagnosis", "staging"),
                    quote = FALSE, noSpaces = TRUE,
                    printToggle = FALSE, showAllLevels = TRUE)
table1 <- cbind(tab1b, tab1_total)
write.csv(table1, "Results/table1.csv")

# How many in each cohort have only one measurement ?
table(df_base$count_long, df_base$cohort)


# Missing values are < 5% in outcome variables. Therefore we can remove these rows without significant bias (only 2 rows are missing the time variable).
df2 <- df1[ rowSums( !is.na(df1[, c("out_fvc", "out_svc", "out_snip",
                                    "out_peak")]) ) ==4 ,  ]


#### What is the raw correlation between outcomes ? ####
M <- cor(df2[, c("out_fvc", "out_svc", "out_peak", "out_snip")])
corrplot(M, method = "number", type = "full")
# What is the raw correlation between outcomes in spinal patients only ?
M_sp <- cor(df2[ df2$diagnosis == "ALS Spinal onset", c("out_fvc", "out_svc", "out_peak", "out_snip")])
corrplot(M_sp, method = "number", type = "full")
M_bulb <- cor(df2[ df2$diagnosis == "ALS Bulbar onset", c("out_fvc", "out_svc", "out_peak", "out_snip")])
corrplot(M_bulb, method = "number", type = "full")



## Bayesian Multiple outcomes longitudinal model
#To determine the agreement across the four outcome metrics, we will use a multiple outcome longitudinal model implemented in Stan via the package brms. In this model, outcomes are considered to be samples from a multivariate normal distribution as a function of independent variables - which may vary per outcome - and with a number of options available to represent the correlation between outcomes. In our modeling we will use the same independent variables for each outcome as they are all metrics of respiratory function. Since our four outcomes are measured contemporaneously in individuals, we will use a random effects longitudinal model. Site of onset will be included as a fixed effect in interaction with time 


#### First perform prior predictive check as per Bayesian workflow ####
fit_priors <- brm(
    mvbind(out_fvc, out_svc, out_snip, out_peak) ~ days_from_baseline * diagnosis +
        gender + cohort + (1 | uin),
    data = df2, chains=8, cores = 8,
    prior = c(set_prior("normal(0, 0.1)", class = "b", resp = c("outfvc", "outsvc","outsnip", "outpeak" )),
              set_prior("normal(3.2, 0.1)", class = "Intercept", resp = "outfvc" ),
              set_prior("normal(3.1, 0.2)", class = "Intercept", resp = "outsvc" ),
              set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
              set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" )
    ),
    sample_prior = "only")


pp1 <- pp_check(fit_priors, resp = "outfvc", nsamples = 200) + ggtitle("FVC") + xlim(-100, 100)
pp2 <- pp_check(fit_priors, resp = "outsvc", nsamples = 200) + ggtitle("SVC") + xlim(-100, 100)
pp3 <- pp_check(fit_priors, resp = "outsnip", nsamples = 200) + ggtitle("SNIP") + xlim(-250, 500)
pp4 <- pp_check(fit_priors, resp = "outpeak", nsamples = 200) + ggtitle("PEAK") + xlim(-1000, 1500)
print(pp1 + pp2 + pp3 + pp4)

tiff("Results/prior_predcheck_models1to3.tiff", width=1000, height=700, res=108)
    print(pp1 + pp2 + pp3 + pp4)
dev.off()



# Prior predictive checks are ok - fit model
fit_mod1 <- brm(
    mvbind(out_fvc, out_svc, out_snip, out_peak) ~ days_from_baseline * diagnosis +
        gender + cohort + (1 | uin),
    data = df2, chains=4, cores = 4, iter = 3500,
    prior = c(set_prior("normal(0, 0.1)", class = "b", resp = c("outfvc", "outsvc","outsnip", "outpeak" )),
              set_prior("normal(3.2, 0.1)", class = "Intercept", resp = "outfvc" ),
              set_prior("normal(3.1, 0.2)", class = "Intercept", resp = "outsvc" ),
              set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
              set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" )
    ),
    sample_prior = "no")
saveRDS(fit_mod1, "Results/brms_fit1_3500.RDS")
#fit_mod1 <- readRDS("Results/brms_fit1_3500.RDS")
summary(fit_mod1)

# trace plot
mcmc_trace(fit_mod1, pars = c("b_outfvc_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outsnip_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outpeak_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outsvc_days_from_baseline:diagnosisALSBulbaronset") )

# Plot prior predictive check as per Bayesian workflow
mod1_pp1 <- pp_check(fit_mod1, resp = "outfvc", nsamples = 200) + ggtitle("FVC")
mod1_pp2 <- pp_check(fit_mod1, resp = "outsvc", nsamples = 200) + ggtitle("SVC")
mod1_pp3 <- pp_check(fit_mod1, resp = "outsnip", nsamples = 200) + ggtitle("SNIP")
mod1_pp4 <- pp_check(fit_mod1, resp = "outpeak", nsamples = 200) + ggtitle("PEAK")

print(mod1_pp1 + mod1_pp2 + mod1_pp3 + mod1_pp4)


# Fit model 2 with more complex correlation structure including time
fit_mod2 <- brm(
    mvbind(out_fvc, out_svc, out_snip, out_peak) ~ days_from_baseline * diagnosis +
        gender + cohort + (days_from_baseline | uin),
    data = df2, chains=4, cores = 4, iter = 3500,
    prior = c(set_prior("normal(0, 0.1)", class = "b", resp = c("outfvc", "outsvc","outsnip", "outpeak" )),
              set_prior("normal(3.2, 0.1)", class = "Intercept", resp = "outfvc" ),
              set_prior("normal(3.1, 0.2)", class = "Intercept", resp = "outsvc" ),
              set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
              set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" )
    ),
    sample_prior = "no")
saveRDS(fit_mod2, "Results/brms_fit2_3500.RDS")
#fit_mod2 <- readRDS("Results/brms_fit2_3500.RDS")
summary(fit_mod2)

mcmc_trace(fit_mod2, pars = c("b_outfvc_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outsnip_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outpeak_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outsvc_days_from_baseline:diagnosisALSBulbaronset") )

# # Perform posterior predictive check as per Bayesian workflow
mod2_pp1 <- pp_check(fit_mod2, resp = "outfvc", nsamples = 200) + ggtitle("FVC")
mod2_pp2 <- pp_check(fit_mod2, resp = "outsvc", nsamples = 200) + ggtitle("SVC")
mod2_pp3 <- pp_check(fit_mod2, resp = "outsnip", nsamples = 200) + ggtitle("SNIP")
mod2_pp4 <- pp_check(fit_mod2, resp = "outpeak", nsamples = 200) + ggtitle("PEAK")
print(mod2_pp1 + mod2_pp2 + mod2_pp3 + mod2_pp4)


# Fit model 3 - extend model 2 to also include group correlation across outcomes
fit_mod3 <- brm(
    mvbind(out_fvc, out_svc, out_snip, out_peak) ~ days_from_baseline * diagnosis +
        gender + cohort + (1 | p | uin) + (0 + days_from_baseline | uin),
    data = df2, chains=4, cores = 4, iter = 3500,
    prior = c(set_prior("normal(0, 0.1)", class = "b", resp = c("outfvc", "outsvc","outsnip", "outpeak" )),
              set_prior("normal(3.2, 0.1)", class = "Intercept", resp = "outfvc" ),
              set_prior("normal(3.1, 0.2)", class = "Intercept", resp = "outsvc" ),
              set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
              set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" )
    ),
    sample_prior = "no")
saveRDS(fit_mod3, "Results/brms_fit3_3500.RDS")
#fit_mod3 <- readRDS("Results/brms_fit3_3500.RDS")


mcmc_trace(fit_mod3, pars = c("b_outfvc_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outsnip_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outpeak_days_from_baseline:diagnosisALSBulbaronset",
                              "b_outsvc_days_from_baseline:diagnosisALSBulbaronset") )

# Perform posterior predictive check as per Bayesian workflow
mod3_pp1 <- pp_check(fit_mod3, resp = "outfvc", nsamples = 200) + ggtitle("FVC")
mod3_pp2 <- pp_check(fit_mod3, resp = "outsvc", nsamples = 200) + ggtitle("SVC")
mod3_pp3 <- pp_check(fit_mod3, resp = "outsnip", nsamples = 200) + ggtitle("SNIP")
mod3_pp4 <- pp_check(fit_mod3, resp = "outpeak", nsamples = 200) + ggtitle("PEAK")
print(mod3_pp1 + mod3_pp2 + mod3_pp3 + mod3_pp4)

tiff("Results/Model3_posterior_predcheck.tiff", width=1000, height=700, res=108)
    print(mod3_pp1 + mod3_pp2 + mod3_pp3 + mod3_pp4)
dev.off()

#extract model fixed effect coefficients and write to file
coefs <- data.frame( fixef(fit_mod3) )
coefs$outcome <- substr( rownames(coefs), 4, 7 )
coefs$variable <- substr( rownames(coefs), 8 , 1000000L)
coefs <- coefs %>% select(outcome, variable, everything())
write.csv(coefs, "Results/model3_unformattedfixed_effects.csv", row.names = FALSE)


# Plot of conditional effects
plotdata <- conditional_effects(fit_mod3, "days_from_baseline:diagnosis")

# make new site onset variable for plotting legend
plotdata$`outfvc.outfvc_days_from_baseline:diagnosis`$site <- factor( ifelse(
    plotdata$`outfvc.outfvc_days_from_baseline:diagnosis`$effect2__ ==
        "ALS Spinal onset", "Spinal", "Bulbar"), levels = c("Spinal", "Bulbar"))
plotdata$`outsvc.outsvc_days_from_baseline:diagnosis`$site <- factor( ifelse(
    plotdata$`outsvc.outsvc_days_from_baseline:diagnosis`$effect2__ ==
        "ALS Spinal onset", "Spinal", "Bulbar"), levels = c("Spinal", "Bulbar"))

plotdata$`outsnip.outsnip_days_from_baseline:diagnosis`$site <- factor( ifelse(
    plotdata$`outsnip.outsnip_days_from_baseline:diagnosis`$effect2__ ==
        "ALS Spinal onset", "Spinal", "Bulbar"), levels = c("Spinal", "Bulbar"))

plotdata$`outpeak.outpeak_days_from_baseline:diagnosis`$site <- factor( ifelse(
    plotdata$`outpeak.outpeak_days_from_baseline:diagnosis`$effect2__ ==
        "ALS Spinal onset", "Spinal", "Bulbar"), levels = c("Spinal", "Bulbar"))

# Make plots
post1 <- ggplot(plotdata$outfvc.outfvc_days_from_baseline,
                aes(x=days_from_baseline, y=estimate__, col=site)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, fill=site), alpha=0.15) +
    coord_cartesian(ylim = c(0, 4)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("Days from baseline") + ylab("FVC (L)")
post2 <- ggplot(plotdata$outsvc.outsvc_days_from_baseline,
                aes(x=days_from_baseline, y=estimate__, col=site)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, fill=site), alpha=0.15) +
    coord_cartesian(ylim = c(0, 4)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("Days from baseline") + ylab("SVC (L)")
post3 <- ggplot(plotdata$outsnip.outsnip_days_from_baseline,
                aes(x=days_from_baseline, y=estimate__, col=site)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, fill=site), alpha=0.15) +
    coord_cartesian(ylim = c(0, 75)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("Days from baseline") + ylab("SNIP (cmH2O )")
post4 <- ggplot(plotdata$outpeak.outpeak_days_from_baseline,
                aes(x=days_from_baseline, y=estimate__, col=site)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, fill=site), alpha=0.15) +
    coord_cartesian(ylim = c(0, 450)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("Days from baseline") + ylab("PEAK (L/min)")

post1 + post2 + post3 + post4 + plot_layout(guides="collect")

tiff("Results/conditional_effects_mod3_figure1.tiff", width=1000, height=700, res=108)
    post1 + post2 + post3 + post4 + plot_layout(guides="collect")
dev.off()

# Extraction variance-covariance matrix
corr_results <- VarCorr(fit_mod3)
write.csv(corr_results, "Results/mod3_varcovar.csv")




##### Compare BRMS models via K-fold cross validation #####
set.seed(45)
kfold_idx2 <- loo::kfold_split_random(K = n_folds1, N = nrow(df2))

# Use parallel computations to fit the folded models
# Note this is computationally demanding and will take several hours for each model
plan(multisession, workers = 10)

k_mod1b <- kfold(fit_mod1, K = n_folds1, folds = kfold_idx2, save_fits = TRUE, chains =1)
saveRDS(k_mod1b, "Results/kfold_mod1b.RDS")
rm(fit_mod1b)
rm(k_mod1b)
gc()

k_mod2b <- kfold(fit_mod2, K = n_folds1, folds = kfold_idx2, save_fits = TRUE, chains =1)
saveRDS(k_mod2b, "Results/kfold_mod2b.RDS")
rm(fit_mod2b)
rm(k_mod2b)
gc()

k_mod3b <- kfold(fit_mod3, K = n_folds1, folds = kfold_idx2, save_fits = TRUE, chains =1)
saveRDS(k_mod3b, "Results/kfold_mod3b.RDS")
rm(fit_mod3b)
rm(k_mod3b)
gc()


# Reload and compare folded models via loo_compare()
k_mod1b <- readRDS("Results/kfold_mod1b.RDS")
k_mod2b <- readRDS("Results/kfold_mod2b.RDS")
k_mod3b <- readRDS("Results/kfold_mod3b.RDS")

loo_compare(k_mod1b, k_mod2b, k_mod3b)

#k_mod3 gives best fit





#### Extract posterior predictions from model at set follow up times ####

# First create new dataset on which to perform predictions
new_data <- df_base[ , c("uin", "cohort", "diagnosis", "gender")]
days_from_baseline <- c(0, 183, 365, 548, 730) # 0, 6, 12, 18, 24 months
new_data <- expand_grid(new_data, days_from_baseline)

# Extract posterior predictions
post_samples <- posterior_predict(fit_mod3, new_data, nsamples = 2000, allow_new_levels = TRUE)
n_samples <- dim(post_samples)[1]

# Create a list of datasets for each sample
post_datasets <- lapply(1:dim(post_samples)[1], function(i) {
    data.frame(new_data, post_samples[ i, ,])
})


# calculate correlation of outcomes overall
cor_all <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][, c("outfvc", "outsvc", "outsnip", "outpeak")]  ))
apply(simplify2array(cor_all), 1:2, mean)


# calculate correlation of outcomes for fixed followup times at 0, 6, 12, 18 and 24 months
idx_time0 <- which( post_datasets[[1]]$days_from_baseline ==0)
cor_time0 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time0, c("outfvc", "outsvc", "outsnip", "outpeak")]  ))

idx_time183 <- which( post_datasets[[1]]$days_from_baseline ==183)
cor_time183 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time183, c("outfvc", "outsvc", "outsnip", "outpeak")]  ))

idx_time365 <- which( post_datasets[[1]]$days_from_baseline ==365)
cor_time365 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time365, c("outfvc", "outsvc", "outsnip", "outpeak")]  ))

idx_time548 <- which( post_datasets[[1]]$days_from_baseline ==548)
cor_time548 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time548, c("outfvc", "outsvc", "outsnip", "outpeak")]  ))

idx_time730 <- which( post_datasets[[1]]$days_from_baseline ==730)
cor_time730 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time730, c("outfvc", "outsvc", "outsnip", "outpeak")]  ))


# Get mean correlations across posterior predictions at each follow-up time
mean_cor_time0 <- apply(simplify2array(cor_time0), 1:2, mean)
mean_cor_time183 <- apply(simplify2array(cor_time183), 1:2, mean)
mean_cor_time365 <- apply(simplify2array(cor_time365), 1:2, mean)
mean_cor_time548 <- apply(simplify2array(cor_time548), 1:2, mean)
mean_cor_time730 <- apply(simplify2array(cor_time730), 1:2, mean)


# Write results to file
write.csv(mean_cor_time0, "Results/model_corr_time0.csv")
write.csv(mean_cor_time183, "Results/model_corr_time183.csv")
write.csv(mean_cor_time365, "Results/model_corr_time365.csv")
write.csv(mean_cor_time548, "Results/model_corr_time548.csv")
write.csv(mean_cor_time730, "Results/model_corr_time730.csv")


#### Extract slopes of fixed effects from fit_mod3 ####

pop_effects <- data.frame( fixef(fit_mod3) )
pop_effects$par <- rownames(pop_effects)
rownames(pop_effects) <- NULL
pop_effects <- pop_effects %>% select(par, everything())

# Remove unwanted pars
wanted_pars <- c("outfvc_Intercept", "outsvc_Intercept", "outsnip_Intercept", 
                 "outpeak_Intercept", "outfvc_days_from_baseline", "outfvc_diagnosisALSBulbaronset", 
                 "outfvc_days_from_baseline:diagnosisALSBulbaronset", 
                 "outsvc_days_from_baseline", "outsvc_diagnosisALSBulbaronset", 
                 "outsvc_days_from_baseline:diagnosisALSBulbaronset", 
                 "outsnip_days_from_baseline", "outsnip_diagnosisALSBulbaronset", 
                 "outsnip_days_from_baseline:diagnosisALSBulbaronset", "outpeak_days_from_baseline", 
                 "outpeak_diagnosisALSBulbaronset",
                 "outpeak_days_from_baseline:diagnosisALSBulbaronset"
)
pop_effects <- pop_effects %>% filter(par %in% c(wanted_pars))
pop_effects[ pop_effects$par == "outfvc_days_from_baseline", ]

# Make a variable called par with reader-friendly parameter names
pop_effects$par <- case_when(pop_effects$par == "outfvc_Intercept" ~ "FVC spinal Intercept",
              pop_effects$par == "outsvc_Intercept" ~ "SVC spinal Intercept",
              pop_effects$par == "outsnip_Intercept" ~ "SNIP spinal Intercept",
              pop_effects$par == "outpeak_Intercept" ~ "PEAK spinal Intercept",
              pop_effects$par == "outfvc_days_from_baseline" ~ "FVC spinal slope",
              pop_effects$par == "outfvc_diagnosisALSBulbaronset" ~ "FVC bulbar Intercept offset",
              pop_effects$par == "outfvc_days_from_baseline:diagnosisALSBulbaronset" ~
                  "FVC bulbar slope offset",
              pop_effects$par == "outsvc_days_from_baseline" ~ "SVC spinal slope",
              pop_effects$par == "outsvc_diagnosisALSBulbaronset" ~ "SVC bulbar Intercept offset",
              pop_effects$par == "outsvc_days_from_baseline:diagnosisALSBulbaronset" ~
                  "SVC bulbar slope offset",
              pop_effects$par == "outsnip_days_from_baseline" ~ "SNIP spinal slope",
              pop_effects$par == "outsnip_diagnosisALSBulbaronset" ~ "SNIP bulbar Intercept offset",
              pop_effects$par == "outsnip_days_from_baseline:diagnosisALSBulbaronset" ~
                  "SNIP bulbar slope offset",
              pop_effects$par == "outpeak_days_from_baseline" ~ "PEAK spinal slope",
              pop_effects$par == "outpeak_diagnosisALSBulbaronset" ~ "PEAK bulbar Intercept offset",
              pop_effects$par == "outpeak_days_from_baseline:diagnosisALSBulbaronset" ~
                  "PEAK bulbar slope offset")

# Filter pop-effects for desired data
FVC_bulb_int <- pop_effects %>% filter(par == "FVC bulbar Intercept offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "FVC spinal Intercept") %>% select(2:5)
FVC_bulb_slope <- pop_effects %>% filter(par == "FVC bulbar slope offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "FVC spinal slope") %>% select(2:5)

SVC_bulb_int <- pop_effects %>% filter(par == "SVC bulbar Intercept offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "SVC spinal Intercept") %>% select(2:5)
SVC_bulb_slope <- pop_effects %>% filter(par == "SVC bulbar slope offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "SVC spinal slope") %>% select(2:5)

SNIP_bulb_int <- pop_effects %>% filter(par == "SNIP bulbar Intercept offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "SNIP spinal Intercept") %>% select(2:5)
SNIP_bulb_slope <- pop_effects %>% filter(par == "SNIP bulbar slope offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "SNIP spinal slope") %>% select(2:5)

PEAK_bulb_int <- pop_effects %>% filter(par == "PEAK bulbar Intercept offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "PEAK spinal Intercept") %>% select(2:5)
PEAK_bulb_slope <- pop_effects %>% filter(par == "PEAK bulbar slope offset") %>% select(2:5) + 
    pop_effects %>% filter(par == "PEAK spinal slope") %>% select(2:5)

FVC_bulb_int$par <- "FVC bulbar Intercept"
FVC_bulb_slope$par <- "FVC bulbar slope"
SVC_bulb_int$par <- "SVC bulbar Intercept"
SVC_bulb_slope$par <- "SVC bulbar slope"
SNIP_bulb_int$par <- "SNIP bulbar Intercept"
SNIP_bulb_slope$par <- "SNIP bulbar slope"
PEAK_bulb_int$par <- "PEAK bulbar Intercept"
PEAK_bulb_slope$par <- "PEAK bulbar slope"

# Bundle extracted efffects together
pop_effects <- bind_rows(pop_effects, FVC_bulb_int, FVC_bulb_slope,
                         SVC_bulb_int, SVC_bulb_slope, SNIP_bulb_int,
                         SNIP_bulb_slope, PEAK_bulb_int, PEAK_bulb_slope)

# Break reader friendly parameter names into multiple columns
metadata <- pop_effects$par %>% str_split(" ", simplify = TRUE) %>% data.frame()
names(metadata) <- c("outcome", "onset", "param", "offsetpar")
pop_effects <- bind_cols(metadata, pop_effects)
pop_effects$par <- NULL

pop_effects <- pop_effects %>% filter( offsetpar!= "offset")
pop_effects$offsetpar <- NULL

# timescale of slopes currently slopes per day - convert to per average month
pop_effects [ pop_effects$param == "slope", 4:7] <- (365.25 *
    pop_effects [ pop_effects$param == "slope", 4:7]) / 12

pop_effects %>% arrange(outcome, param, onset)

# Write to file
write.csv(pop_effects, "Results/model3_Interceptsandslopes.csv")


#### Refit model3 using percent predicted outcomes ####

# Recheck prior predictive check with new prior values to reflect changed outcome scale
fit_priors2 <- brm(
    mvbind(fvc_max_percent_pred, svc_percent_pred, out_snip, out_peak) ~ days_from_baseline * diagnosis +
        gender + cohort + (1 | uin),
    data = df2, chains=8, cores = 8,
    prior = c(set_prior("normal(0, 0.1)", class = "b", resp = c("fvcmaxpercentpred", "svcpercentpred","outsnip", "outpeak" )),
              set_prior("normal(75, 5)", class = "Intercept", resp = "fvcmaxpercentpred" ),
              set_prior("normal(75, 5)", class = "Intercept", resp = "svcpercentpred" ),
              set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
              set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" )
    ),
    sample_prior = "only")

pp1 <- pp_check(fit_priors2, resp = "fvcmaxpercentpred", nsamples = 200) + ggtitle("FVC") + xlim(-100, 200)
pp2 <- pp_check(fit_priors2, resp = "svcpercentpred", nsamples = 200) + ggtitle("SVC") + xlim(-100, 200)
pp3 <- pp_check(fit_priors2, resp = "outsnip", nsamples = 200) + ggtitle("SNIP") + xlim(-250, 500)
pp4 <- pp_check(fit_priors2, resp = "outpeak", nsamples = 200) + ggtitle("PEAK") + xlim(-1000, 1500)
print(pp1 + pp2 + pp3 + pp4)

tiff("Results/prior_predcheck_percentscale.tiff", width=1000, height=700, res=108)
    print(pp1 + pp2 + pp3 + pp4)
dev.off()



# Fit model 3b - same explanatory variables as model 3, but with FVC & SVC outcomes on precent predicted scale
fit_mod3b <- brm(
    mvbind(fvc_max_percent_pred, svc_percent_pred, out_snip, out_peak) ~ days_from_baseline * diagnosis +
        gender + cohort + (1 | p | uin) + (0 + days_from_baseline | uin),
    data = df2, chains=4, cores = 4, iter = 3500,
    prior = c(set_prior("normal(0, 0.1)", class = "b", resp = c("fvcmaxpercentpred", "svcpercentpred","outsnip", "outpeak" )),
              set_prior("normal(75, 5)", class = "Intercept", resp = "fvcmaxpercentpred" ),
              set_prior("normal(75, 5)", class = "Intercept", resp = "svcpercentpred" ),
              set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
              set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" )
    ),
    sample_prior = "no")
saveRDS(fit_mod3b, "Results/brms_fit3b_3500.RDS")
#fit_mod3b <- readRDS( "Results/brms_fit3b_3500.RDS")

# Check traceplot
mcmc_trace(fit_mod3b, pars = c("b_fvcmaxpercentpred_days_from_baseline:diagnosisALSBulbaronset",
                               "b_outsnip_days_from_baseline:diagnosisALSBulbaronset",
                               "b_outpeak_days_from_baseline:diagnosisALSBulbaronset",
                               "b_svcpercentpred_days_from_baseline:diagnosisALSBulbaronset") )

# Perform posterior predictive check as per Bayesian workflow
mod3b_pp1 <- pp_check(fit_mod3b, resp = "fvcmaxpercentpred", nsamples = 200) + ggtitle("FVC")
mod3b_pp2 <- pp_check(fit_mod3b, resp = "svcpercentpred", nsamples = 200) + ggtitle("SVC")
mod3b_pp3 <- pp_check(fit_mod3b, resp = "outsnip", nsamples = 200) + ggtitle("SNIP")
mod3b_pp4 <- pp_check(fit_mod3b, resp = "outpeak", nsamples = 200) + ggtitle("PEAK")
print(mod3b_pp1 + mod3b_pp2 + mod3b_pp3 + mod3b_pp4)

# Make plot of conditional effects
plotdata3b <- conditional_effects(fit_mod3b, "days_from_baseline:diagnosis")

post1_3b <- ggplot(plotdata3b$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:diagnosis` ,
                   aes(x=days_from_baseline, y=estimate__, col=effect2__)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA), alpha=0.15) +
    coord_cartesian(ylim = c(0, 100)) + theme_minimal() + labs(color='Site of Onset') +
    xlab("Days from baseline") + ylab("FVC % of max predicted")
post2_3b <- ggplot(plotdata3b$`svcpercentpred.svcpercentpred_days_from_baseline:diagnosis`,
                   aes(x=days_from_baseline, y=estimate__, col=effect2__)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA), alpha=0.15) +
    coord_cartesian(ylim = c(0, 100)) + theme_minimal() + labs(color='Site of Onset') +
    xlab("Days from baseline") + ylab("SVC % of max predicted")
post3_3b <- ggplot(plotdata3b$outsnip.outsnip_days_from_baseline,
                   aes(x=days_from_baseline, y=estimate__, col=effect2__)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA), alpha=0.15) +
    coord_cartesian(ylim = c(0, 75)) + theme_minimal() + labs(color='Site of Onset') +
    xlab("Days from baseline") + ylab("SNIP")
post4_3b <- ggplot(plotdata3b$outpeak.outpeak_days_from_baseline,
                   aes(x=days_from_baseline, y=estimate__, col=effect2__)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA), alpha=0.15) +
    coord_cartesian(ylim = c(0, 450)) + theme_minimal() + labs(color='Site of Onset') +
    xlab("Days from baseline") + ylab("PEAK")

post1_3b + post2_3b + post3_3b + post4_3b + plot_layout(guides="collect")


# Extract slopes of fixed effects from fit_mod3b
pop_effects3b <- data.frame( fixef(fit_mod3b) )
pop_effects3b$par <- rownames(pop_effects3b)
rownames(pop_effects3b) <- NULL
pop_effects3b <- pop_effects3b %>% select(par, everything())

# Remove unwanted pars
wanted_pars <- c("fvcmaxpercentpred_Intercept", "svcpercentpred_Intercept", "outsnip_Intercept", 
                 "outpeak_Intercept", "fvcmaxpercentpred_days_from_baseline", "fvcmaxpercentpred_diagnosisALSBulbaronset", 
                 "fvcmaxpercentpred_days_from_baseline:diagnosisALSBulbaronset", 
                 "svcpercentpred_days_from_baseline", "svcpercentpred_diagnosisALSBulbaronset", 
                 "svcpercentpred_days_from_baseline:diagnosisALSBulbaronset", 
                 "outsnip_days_from_baseline", "outsnip_diagnosisALSBulbaronset", 
                 "outsnip_days_from_baseline:diagnosisALSBulbaronset", "outpeak_days_from_baseline", 
                 "outpeak_diagnosisALSBulbaronset",
                 "outpeak_days_from_baseline:diagnosisALSBulbaronset"
)
pop_effects3b <- pop_effects3b %>% filter(par %in% c(wanted_pars))
pop_effects3b[ pop_effects3b$par == "fvcmaxpercentpred_days_from_baseline", ]

# Make a variable called par with reader-friendly parameter names
pop_effects3b$par <- case_when(pop_effects3b$par == "fvcmaxpercentpred_Intercept" ~ "FVC% spinal Intercept",
                               pop_effects3b$par == "svcpercentpred_Intercept" ~ "SVC% spinal Intercept",
                               pop_effects3b$par == "outsnip_Intercept" ~ "SNIP spinal Intercept",
                               pop_effects3b$par == "outpeak_Intercept" ~ "PEAK spinal Intercept",
                               pop_effects3b$par == "fvcmaxpercentpred_days_from_baseline" ~ "FVC% spinal slope",
                               pop_effects3b$par == "fvcmaxpercentpred_diagnosisALSBulbaronset" ~ "FVC% bulbar Intercept offset",
                               pop_effects3b$par == "fvcmaxpercentpred_days_from_baseline:diagnosisALSBulbaronset" ~
                                   "FVC% bulbar slope offset",
                               pop_effects3b$par == "svcpercentpred_days_from_baseline" ~ "SVC% spinal slope",
                               pop_effects3b$par == "svcpercentpred_diagnosisALSBulbaronset" ~ "SVC% bulbar Intercept offset",
                               pop_effects3b$par == "svcpercentpred_days_from_baseline:diagnosisALSBulbaronset" ~
                                   "SVC% bulbar slope offset",
                               pop_effects3b$par == "outsnip_days_from_baseline" ~ "SNIP spinal slope",
                               pop_effects3b$par == "outsnip_diagnosisALSBulbaronset" ~ "SNIP bulbar Intercept offset",
                               pop_effects3b$par == "outsnip_days_from_baseline:diagnosisALSBulbaronset" ~
                                   "SNIP bulbar slope offset",
                               pop_effects3b$par == "outpeak_days_from_baseline" ~ "PEAK spinal slope",
                               pop_effects3b$par == "outpeak_diagnosisALSBulbaronset" ~ "PEAK bulbar Intercept offset",
                               pop_effects3b$par == "outpeak_days_from_baseline:diagnosisALSBulbaronset" ~
                                   "PEAK bulbar slope offset")

# Filter pop-effects for desired data
FVC_bulb_int <- pop_effects3b %>% filter(par == "FVC% bulbar Intercept offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "FVC% spinal Intercept") %>% select(2:5)
FVC_bulb_slope <- pop_effects3b %>% filter(par == "FVC% bulbar slope offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "FVC% spinal slope") %>% select(2:5)

SVC_bulb_int <- pop_effects3b %>% filter(par == "SVC% bulbar Intercept offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "SVC% spinal Intercept") %>% select(2:5)
SVC_bulb_slope <- pop_effects3b %>% filter(par == "SVC% bulbar slope offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "SVC% spinal slope") %>% select(2:5)

SNIP_bulb_int <- pop_effects3b %>% filter(par == "SNIP bulbar Intercept offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "SNIP spinal Intercept") %>% select(2:5)
SNIP_bulb_slope <- pop_effects3b %>% filter(par == "SNIP bulbar slope offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "SNIP spinal slope") %>% select(2:5)

PEAK_bulb_int <- pop_effects3b %>% filter(par == "PEAK bulbar Intercept offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "PEAK spinal Intercept") %>% select(2:5)
PEAK_bulb_slope <- pop_effects3b %>% filter(par == "PEAK bulbar slope offset") %>% select(2:5) + 
    pop_effects3b %>% filter(par == "PEAK spinal slope") %>% select(2:5)

FVC_bulb_int$par <- "FVC% bulbar Intercept"
FVC_bulb_slope$par <- "FVC% bulbar slope"
SVC_bulb_int$par <- "SVC% bulbar Intercept"
SVC_bulb_slope$par <- "SVC% bulbar slope"
SNIP_bulb_int$par <- "SNIP bulbar Intercept"
SNIP_bulb_slope$par <- "SNIP bulbar slope"
PEAK_bulb_int$par <- "PEAK bulbar Intercept"
PEAK_bulb_slope$par <- "PEAK bulbar slope"

# Bundle extracted efffects together
pop_effects3b <- bind_rows(pop_effects3b, FVC_bulb_int, FVC_bulb_slope,
                           SVC_bulb_int, SVC_bulb_slope, SNIP_bulb_int,
                           SNIP_bulb_slope, PEAK_bulb_int, PEAK_bulb_slope)

# Break reader friendly parameter names into multiple columns
metadata <- pop_effects3b$par %>% str_split(" ", simplify = TRUE) %>% data.frame()
names(metadata) <- c("outcome", "onset", "param", "offsetpar")

pop_effects3b <- bind_cols(metadata, pop_effects3b)
pop_effects3b$par <- NULL

pop_effects3b <- pop_effects3b %>% filter( offsetpar!= "offset")
pop_effects3b$offsetpar <- NULL

# Slopes currently per day - convert to per average month
pop_effects3b [ pop_effects3b$param == "slope", 4:7] <- (365.25 *
                                                             pop_effects3b [ pop_effects3b$param == "slope", 4:7]) / 12

pop_effects3b %>% arrange(outcome, param, onset)

# Write results to file
write.csv(pop_effects3b, "Results/percentfit_Interceptsandslopes.csv")
