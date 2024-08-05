### Drivers of FKBP5 DNA methylation: source functions ###
# Author: Vera N. Karlbauer
# Contact: vera_karlbauer@psych.mpg.de
# Date created: 25-07-20204
# Purpose: Source functions for differential DNA methylation analysis with predictor(s) of interest, used in scripts 05-08.

### Function 1: Loop simple linear model over 49 GpGs & return results

#' loop_lm_simple
#'
#' @param data a dataframe of methylation data and covariates with subject IDs as rows. Important: methylation data takes up first 49 columns, otherwise this will not work.
#' @param predictor_name character value specifiying column name in data for the main predictor of interest for which p values are returned
#' @param object_name character value specifiying name for the final dataframe with results
#'
#' @return dataframe with adjusted R, partial R, coefficients, p-values, fdr-adjusted p-values, and significance as "yes"/"no" per CpG site
#' @export
#'
#' @examples
#' # Differential methylation for genotype on dataset all_t0
#' loop_lm_simple(data = all_t0, predictor_name = "genotype", object_name = "diffmeth_geno_t0")

loop_lm_simple <- function(data, predictor_name){
  # set up results object with CpG names
  resultsobj <- cbind(colnames(data[1:49]),colnames(data[1:49]),colnames(data[1:49]),colnames(data[1:49]),colnames(data[1:49]))
  # loop linear models over all GpGs
  for(i in 1:49){
    # linear model with pre-specified covariates
    # age, sex
    fit = lm(formula = data[,i] ~ data$age + data$Sex 
             # genetic components
             + data$C1 + data$C2 + data$C3 
             # cell type counts
             + data$CD14 + data$CD34 + data$Buccal 
             # prenatal smoke & alcohol exposure
             + data$smoke_exposure_score + data$alcohol_exposure_score 
             # predictor of interest
             + data[, predictor_name])
    # extract summary measures
    s = summary(fit)
    rsqpar = suppressWarnings(rsq.partial(fit, adj = TRUE))
    resultsobj[i,2] = s$adj.r.squared
    resultsobj[i,3] = s$coefficients["data[, predictor_name]", "Pr(>|t|)"]
    resultsobj[i,4] = s$coefficients["data[, predictor_name]", "Estimate"]
    resultsobj[i,5] = rsqpar$partial.rsq[11]
  }
  # results object as data frame
  resultsobj <- as.data.frame(resultsobj)
  # rename columns in results onjext
  colnames(resultsobj) = c("CpG","adj.R","p", "coef", "partial_adj.R")
  # set to numeric
  resultsobj <- resultsobj %>%
    mutate_at(c("adj.R", "p", "coef", "partial_adj.R"), as.numeric)
  # add column calling uncorrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_uncorrected = case_when(
                         p >= 0.05 ~ "no",
                         p < 0.05 ~ "yes")
  )
  # correct p-values
  resultsobj$p_adjust <- as.numeric(p.adjust(p = resultsobj$p, method = "fdr"))
  # add columnn calling corrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_corrected = case_when(
                         p_adjust >= 0.05 ~ "no",
                         p_adjust < 0.05 ~ "yes"))
                       # relocate columns
   resultsobj <- relocate(resultsobj,partial_adj.R, .after = adj.R)
  # return results object
  return(resultsobj)
  # print results object to console
  print(resultsobj)
}

##@ Function 2: loop_lm_additive

#' Loop linear model with additive interaction over 49 GpGs & return results
#'
#' @param data a dataframe of methylation data and covariates with subject IDs as rows. Important: methylation data takes up first 49 columns, otherwise this will not work.
#' @param predictor_1 character value specifiying column name in data for the first predictor for the additive interaction of interest for which p values are returned
#' @param predictor_2 character value specifiying column name in data for the first predictor for the additive interaction of interest for which p values are returned
#' @param object_name character value specifiying name for the final dataframe with results
#'
#' @return dataframe with adjusted R, partial R for predictors, coefficients, p-values, fdr-adjusted p-values, and significance as "yes"/"no" per CpG site
#' @export
#'
#' @examples
#' # Additive interaction of genotype and glucocorticoid exposure on methylation on dataset all_t0
#' loop_lm_additive(data = all_t0, predictor_1 = "genotype", predictor_2 = "gc_exposure", object_name = "diffmeth_geno_gc_add_t0")

loop_lm_additive <- function(data, predictor_1, predictor_2){
  # set up results object with CpG names
  resultsobj <- cbind(colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]))
  # loop linear models over all GpGs
  for (i in 1:49){
    # setup null model with all covariates and predictor 1
    # age, sex
    null = lm(formula = data[,i] ~ data$age + data$Sex 
              # genetic components
              + data$C1 + data$C2 + data$C3 
              # cell type counts
              + data$CD14 + data$CD34 + data$Buccal 
              # prenatal smoke & alcohol exposure
              + data$smoke_exposure_score + data$alcohol_exposure_score 
              # predictor of interest
              + data[, predictor_1])
    # set up additive model with all covariates and predictor 1 & 2
    # age, sex
    fit = lm(formula = data[,i] ~ data$age + data$Sex 
             # genetic components
             + data$C1 + data$C2 + data$C3 
             # cell type counts
             + data$CD14 + data$CD34 + data$Buccal 
             # prenatal smoke & alcohol exposure
             + data$smoke_exposure_score + data$alcohol_exposure_score 
             # predictor of interest
             + data[, predictor_1] 
             + data[, predictor_2])
    comparison_additive = anova(null,fit)
    s = summary(fit)
    rsqpar = suppressWarnings(rsq.partial(fit, adj = TRUE))
    pval = comparison_additive[2,6]
    resultsobj[i,2] = s$adj.r.squared
    resultsobj[i,3] = pval
    resultsobj[i,4] = fit$coefficients["data[, predictor_1]"]
    resultsobj[i,5] = fit$coefficients["data[, predictor_2]"]
    resultsobj[i,6] = rsqpar$partial.rsq[11]
    resultsobj[i,7] = rsqpar$partial.rsq[12]
  }
  # results object as data frame
  resultsobj <- as.data.frame(resultsobj)
  # rename columns in results onjext
  colnames(resultsobj) = c("CpG","adj.R","p", "coef_1", "coef_2", "partial_adj.R_1", "partial_adj.R_2")
  # set to numeric
  resultsobj <- resultsobj %>%
    mutate_at(c("adj.R", "p", "coef_1", "coef_2", "partial_adj.R_1", "partial_adj.R_2"), as.numeric)
  # add column calling uncorrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_uncorrected = case_when(
                         p >= 0.05 ~ "no",
                         p < 0.05 ~ "yes")
  )
  # correct p-values
  resultsobj$p_adjust <- as.numeric(p.adjust(p = resultsobj$p, method = "fdr"))
  # add columnn calling corrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_corrected = case_when(
                         p_adjust >= 0.05 ~ "no",
                         p_adjust < 0.05 ~ "yes"))
  # relocate columns
  resultsobj <- resultsobj %>%
    relocate(partial_adj.R_1, .after = adj.R) %>%
    relocate(partial_adj.R_2, .after = partial_adj.R_1)
  # rename coefficient columns after predictor names
  # set colnames
  coefnames <- c("CpG", "adj.R", 
                 paste("partial_adj.R", predictor_1, sep = "_"), 
                 paste("partial_adj.R", predictor_2, sep = "_"), 
                 "p", 
                 paste("coef", predictor_1, sep = "_"), 
                 paste("coef", predictor_2, sep = "_"), 
                 "sig_uncorrected", "p_adjust", "sig_corrected")
  # rename
  colnames(resultsobj) <- coefnames
  # name results object and write to working directory
  return(resultsobj)
  # print results object to console
  print(resultsobj)
}

### Function 3: loop_lm_multiplicative

#' Loop linear model with multiplicative interaction over 49 GpGs & return results
#'
#' @param data a dataframe of methylation data and covariates with subject IDs as rows. Important: methylation data takes up first 49 columns, otherwise this will not work.
#' @param predictor_1 character value specifiying column name in data for the first predictor for the multiplicative interaction of interest for which p values are returned
#' @param predictor_2 character value specifiying column name in data for the first predictor for the multiplicative interaction of interest for which p values are returned
#' @param object_name character value specifiying name for the final dataframe with results
#'
#' @return dataframe with adjusted R, coefficients, p-values, fdr-adjusted p-values, and significance as "yes"/"no" per CpG site
#' @export
#'
#' @examples
#' # Multiplicative interaction of genotype and glucocorticoid exposure on methylation on dataset all_t0
#' loop_lm_additive(data = all_t0, predictor_1 = "genotype", predictor_2 = "gc_exposure", object_name = "diffmeth_geno_gc_mult_t0")


loop_lm_multiplicative <- function(data, predictor_1, predictor_2){
  # set up results object with CpG names
  resultsobj <- cbind(colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]), colnames(data[1:49]))
  # loop linear models over all GpGs
  for (i in 1:49){
    # setup null model with all covariates, predictor 1 + predictor 2
    # age, sex
    null = lm(formula = data[,i] ~ data$age + data$Sex 
              # genetic components
              + data$C1 + data$C2 + data$C3 
              # cell type counts
              + data$CD14 + data$CD34 + data$Buccal 
              # prenatal smoke & alcohol exposure
              + data$smoke_exposure_score + data$alcohol_exposure_score 
              # predictors of interest
              + data[, predictor_1] 
              + data[, predictor_2])
    # set up additive model with all covariates and predictor 1 & 2
    # age, sex
    fit = lm(formula = data[,i] ~ data$age + data$Sex 
             # genetic components
             + data$C1 + data$C2 + data$C3 
             # cell type counts
             + data$CD14 + data$CD34 + data$Buccal 
             # prenatal smoke & alcohol exposure
             + data$smoke_exposure_score + data$alcohol_exposure_score 
             # interaction between of predictors of interest
             + data[, predictor_1] * data[, predictor_2])
    comparison_additive = anova(null,fit)
    s = summary(fit)
    rsqpar = suppressWarnings(rsq.partial(fit, adj = TRUE))
    pval = comparison_additive[2,6]
    resultsobj[i,2] = s$adj.r.squared
    resultsobj[i,3] = pval
    resultsobj[i,4] = fit$coefficients["data[, predictor_1]"]
    resultsobj[i,5] = fit$coefficients["data[, predictor_2]"]
    resultsobj[i,6] = fit$coefficients["data[, predictor_1]:data[, predictor_2]"]
    resultsobj[i,7] = rsqpar$partial.rsq[11]
    resultsobj[i,8] = rsqpar$partial.rsq[12]
    resultsobj[i,9] = rsqpar$partial.rsq[13]
  }
  # results object as data frame
  resultsobj <- as.data.frame(resultsobj)
  # rename columns in results onjext
  colnames(resultsobj) = c("CpG","adj.R","p", "coef_1", "coef_2", "coef_mult", "partial_adj.R_1", "partial_adj.R_2", "partial_adj.R_mult")
  # set to numeric
  resultsobj <- resultsobj %>%
    mutate_at(c("adj.R", "p", "coef_1", "coef_2", "coef_mult", "partial_adj.R_1", "partial_adj.R_2", "partial_adj.R_mult"), as.numeric)
  # add column calling uncorrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_uncorrected = case_when(
                         p >= 0.05 ~ "no",
                         p < 0.05 ~ "yes")
  )
  # correct p-values
  resultsobj$p_adjust <- as.numeric(p.adjust(p = resultsobj$p, method = "fdr"))
  # add columnn calling corrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_corrected = case_when(
                         p_adjust >= 0.05 ~ "no",
                         p_adjust < 0.05 ~ "yes"))
  # relocate columns
  resultsobj <- resultsobj %>%
    relocate(partial_adj.R_1, .after = adj.R) %>%
    relocate(partial_adj.R_2, .after = partial_adj.R_1) %>%
    relocate(partial_adj.R_mult, .after = partial_adj.R_2)
  # rename coefficient columns after predictor names
  # set colnames
  coefnames <- c("CpG", "adj.R", 
                 paste("partial_adj.R", predictor_1, sep = "_"), 
                 paste("partial_adj.R", predictor_2, sep = "_"), 
                 "partial_adj.R_mult", "p", 
                 paste("coef", predictor_1, sep = "_"), 
                 paste("coef", predictor_2, sep = "_"), 
                "coef_mult", "sig_uncorrected", "p_adjust", "sig_corrected")
  # rename
  colnames(resultsobj) <- coefnames
  # name results object and write to working directory
  return(resultsobj)
  # print results object to console
  print(resultsobj)
}

### Function 4: loop_lm_time

#' Loop linear model with time effect over 49 GpGs & return results
#'
#' @param data a dataframe of methylation data and covariates with subject IDs as rows. Needs to be coded in long format, i.e. every timepoint has its own row and needs to be specified with the variable "Time. Important: methylation data takes up first 49 columns, otherwise this will not work.
#' @param predictor_name character value specifiying column name in data for the main predictor of interest for which p values are returned
#' @param object_name character value specifiying name for the final dataframe with results
#'
#' @return dataframe with R squared, coefficients, p-values, fdr-adjusted p-values, and significance as "yes"/"no" per CpG site
#' @export
#'
#' @examples 
#' Genotype effects over time
#' loop_lm_time(data = all_time, predictor_name = "genotype, object_name = "diffmeth_geno_time)

loop_lm_time <- function(data, predictor_name){
  # set up results object with CpG names
  resultsobj <- cbind(colnames(data[1:49]),colnames(data[1:49]),colnames(data[1:49]),colnames(data[1:49]))
  for (i in 1:49){
    # subject as factor and back to numeric again
    x = as.numeric(as.factor(data$Subject))
    # setup null model with all covariates, time and subject effect
    # age, sex
    null = lmer(formula = data[,i] ~ data$age + data$Sex 
                # genetic components
                + data$C1 + data$C2 + data$C3 
                # cell type counts
                + data$CD14 + data$CD34 + data$Buccal 
                # prenatal smoke & alcohol exposure
                + data$smoke_exposure_score + data$alcohol_exposure_score 
                # time effect
                + data$Time
                # subject effect
                + (1|x), data = data, REML = F)
    # set up model with all covariates, time and subject, and predictor of interest
    # age, sex
    fit = lmer(formula = data[,i] ~ data$age + data$Sex 
               # genetic components
               + data$C1 + data$C2 + data$C3 
               # cell type counts
               + data$CD14 + data$CD34 + data$Buccal 
               # prenatal smoke & alcohol exposure
               + data$smoke_exposure_score + data$alcohol_exposure_score 
               # time effect
               + data$Time
               # predictor of interest
               + data[, predictor_name] 
               # subject effect
               + (1|x), data = data, REML = F)
    comparison_additive = anova(null,fit)
    s = summary(fit)
    pval = comparison_additive[2,8]
    resultsobj[i,2] = s$AICtab["AIC"]
    resultsobj[i,3] = pval
    resultsobj[i,4] = s$coefficients["data[, predictor_name]", "Estimate"]
  }
  # results object as data frame
  resultsobj <- as.data.frame(resultsobj)
  # rename columns in results object
  colnames(resultsobj) = c("CpG","AIC","p", "coef")
  # set to numeric
  resultsobj <- resultsobj %>%
    mutate_at(c("AIC", "p", "coef"), as.numeric)
  # add column calling uncorrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_uncorrected = case_when(
                         p >= 0.05 ~ "no",
                         p < 0.05 ~ "yes")
  )
  # correct p-values
  resultsobj$p_adjust <- as.numeric(p.adjust(p = resultsobj$p, method = "fdr"))
  # add columnn calling corrected significance
  resultsobj <- mutate(resultsobj, 
                       sig_corrected = case_when(
                         p_adjust >= 0.05 ~ "no",
                         p_adjust < 0.05 ~ "yes")
  )
  # name results object and write to working directory
  return(resultsobj)
  # print results object to console
  print(resultsobj)
}
