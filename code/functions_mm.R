##########################################
####  Functions for mixed modelling
##########################################
#### | Project name: Bioclimatic envelopes
#### | Creator: Mirza Cengic
#### | Contact: mirzaceng@gmail.com
##########################################

####
# Fit mixed models
fit_mixed2 <- function(data, full = TRUE, basic = FALSE)
{
  cat("Simple glmm's", "\n")
  baseline <- lmer(value ~ 1 + (1 | binomial), data = data, REML = FALSE)
  pred_model <- lmer(value ~ 1 + (1 | binomial) + predictor_set, data = data, REML = FALSE)
  algo_model <- lmer(value ~ 1 + (1|binomial) + predictor_set + algorithm, data = data, REML = FALSE)
  pa_model <- lmer(value ~ 1 + (1|binomial) + predictor_set + algorithm + pseudoabs_set, data = data, REML = FALSE)
  niche_model <- lmer(value ~ 1 + (1|binomial) + predictor_set + algorithm + pseudoabs_set + niche_overlap_value, data = data, REML = FALSE)
  cat("With interactions", "\n")
  pred_algoM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm, data = data, REML = FALSE)
  pred_psM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * pseudoabs_set, data = data, REML = FALSE)
  algo_psM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * predictor_set + algorithm * predictor_set, data = data, REML = FALSE)
  
  algo_niche <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * predictor_set + algorithm * predictor_set + algorithm * niche_overlap_value, data = data, REML = FALSE)
  pred_algo_psM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * pseudoabs_set + algorithm * pseudoabs_set + predictor_set * algorithm * predictor_set, data = data, REML = FALSE)
  
  m <- anova(baseline, pred_model, algo_model, pa_model, niche_model, pred_algoM, pred_psM, algo_psM, algo_niche, pred_algo_psM)
  
 
  if (basic == TRUE)
  {
    m1 <- suppressWarnings(broom::tidy(m))
    best_model_name <- m1[which.min(m1$AIC), 1]
    best_mod <- get(as.character(best_model_name))
    return(best_mod)
  }
  
  if (full == FALSE)
  {
    return(broom::tidy(m))
  }
  if (full == TRUE)
  {
    m1 <- broom::tidy(m)
    best_model_name <- m1[which.min(m1$AIC), 1]
    # best_model_name <- pull(m1[which.min(m1$AIC), 1])
    print(best_model_name)
    str(best_model_name)
    r2 <- r2beta(get(as.character(best_model_name)))
    return(r2)
  }
}

get_significance_vals <- function(dat, grp)
{
  my_mod <- dat %>% 
    mutate(
      algorithm = if_else(algorithm == "EMmeanByTSS", "EM", algorithm)
    ) %>% 
    filter(evaluation_set == "Eval") %>% 
    group_by(group, predictor_set, pseudoabs_set, algorithm, evaluation_mode, binomial) %>% 
    distinct(binomial, .keep_all = TRUE) %>% 
    ungroup() %>% 
    filter(group == grp) %>% 
    fit_mixed2(basic = TRUE)
  
  return(my_mod)
}
#### fit_mixed

fit_mixed <- function(data, full = TRUE)
{
  
  baseline <- lmer(value ~ 1 + (1 | binomial), data = data, REML = FALSE)
  pred_model <- lmer(value ~ 1 + (1 | binomial) + predictor_set, data = data, REML = FALSE)
  algo_model <- lmer(value ~ 1 + (1|binomial) + predictor_set + algorithm, data = data, REML = FALSE)
  pa_model <- lmer(value ~ 1 + (1|binomial) + predictor_set + algorithm + pseudoabs_set, data = data, REML = FALSE)
  niche_model <- lmer(value ~ 1 + (1|binomial) + predictor_set + algorithm + pseudoabs_set + niche_overlap_value, data = data, REML = FALSE)
  cat("Simple glmm's", "\n")
  pred_algoM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm, data = data, REML = FALSE)
  pred_psM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * pseudoabs_set, data = data, REML = FALSE)
  algo_psM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * predictor_set + algorithm * predictor_set, data = data, REML = FALSE)
  cat("With interactions", "\n")
  algo_niche <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * predictor_set + algorithm * predictor_set + algorithm * niche_overlap_value, data = data, REML = FALSE)
  pred_algo_psM <- lmer(value ~ (1|binomial) + predictor_set + algorithm + pseudoabs_set + predictor_set * algorithm + predictor_set * pseudoabs_set + algorithm * pseudoabs_set + predictor_set * algorithm * predictor_set, data = data, REML = FALSE)
  
  m <- anova(baseline, pred_model, algo_model, pa_model, niche_model, pred_algoM, pred_psM, algo_psM, algo_niche, pred_algo_psM)
  if (full == FALSE)
  {
    return(broom::tidy(m))
  }
  if (full == TRUE)
  {
    m1 <- suppressWarnings(broom::tidy(m))
    best_model_name <- m1[which.min(m1$AIC), 1]
    # best_model_name <- pull(m1[which.min(m1$AIC), 1])
    print(best_model_name)
    str(best_model_name)
    r2 <- r2beta(get(as.character(best_model_name)))
    return(r2)
  }
}


# Function to process results of mixed models
process_mm_results <- function(data)
{
  data %>% 
    filter(Effect != "Model") %>% 
    mutate(
      Effect = as.character(Effect),
      Effect_grouped = case_when(
        Effect == "algorithmGLM" ~ "Algo",
        Effect == "algorithmGAM" ~ "Algo",
        Effect == "algorithmGBM" ~ "Algo",
        Effect == "algorithmMARS" ~ "Algo",
        Effect == "algorithmMAXENT.Phillips" ~ "Algo",
        Effect == "algorithmMAXENT.Background" ~ "Algo",
        Effect == "algorithmRF" ~ "Algo",
        Effect == "algorithmEnsemble" ~ "Algo",
        Effect == "pseudoabs_setPA=1k" ~ "PseudoAbsence",
        Effect == "pseudoabs_setPA=10k" ~ "PseudoAbsence",
        Effect == "pseudoabs_setPA=P" ~ "PseudoAbsence",
        Effect == "predictor_set2var" ~ "VarSet",
        Effect == "predictor_set4var" ~ "VarSet",
        Effect == "predictor_set10var" ~ "VarSet",
        Effect == "predictor_setallvar" ~ "VarSet",
        Effect == "niche_overlap_value" ~ "Niche",
        Effect == "predictor_setallvar:pseudoabs_setPA=10k" ~ "VarSet*PseudoAbsence",
        Effect == "algorithmGBM:pseudoabs_setPA=10k" ~ "VarSet*Algo",
        Effect == "algorithmRF:niche_overlap_value" ~ "Algo*Niche",
        Effect == "predictor_set4var:algorithmRF" ~ "VarSet*Algo",
        Effect == "predictor_set2var:algorithmRF" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmGLM" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmRF" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmGAM" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmGLM" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmMARS" ~ "VarSet*Algo",
        Effect == "predictor_set2var:algorithmGBM" ~ "VarSet*Algo",
        Effect == "algorithmGLM:niche_overlap_value" ~ "Algo*Niche",
        Effect == "predictor_set2var:algorithmGLM" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmMARS" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmGBM" ~ "VarSet*Algo",
        Effect == "algorithmMARS:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmGBM:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmGAM:niche_overlap_value" ~ "Algo*Niche",
        Effect == "predictor_set2var:algorithmMARS" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmGBM" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmGAM" ~ "VarSet*Algo",
        Effect == "predictor_set2var:pseudoabs_setPA=10k" ~ "VarSet*PseudoAbsence",
        Effect == "algorithmRF:pseudoabs_setPA=10k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:pseudoabs_setPA=10k" ~ "VarSet*PseudoAbsence",
        Effect == "algorithmRF:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:pseudoabs_setPA=P" ~ "VarSet*PseudoAbsence",
        Effect == "algorithmGLM:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmMARS:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmMARS:pseudoabs_setPA=10k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmGLM:pseudoabs_setPA=10k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:pseudoabs_setPA=P" ~ "VarSet*PseudoAbsence",
        Effect == "predictor_set2var:pseudoabs_setPA=P" ~ "VarSet*PseudoAbsence",
        Effect == "predictor_setallvar:pseudoabs_setPA=P" ~ "VarSet*PseudoAbsence",
        Effect == "algorithmGAM:pseudoabs_setPA=10k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmGBM:pseudoabs_setPA=10k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmGBM:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set2var:algorithmGAM" ~ "VarSet*Algo",
        Effect == "algorithmRF:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmRF:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:pseudoabs_setPA=1k" ~ "VarSet*PseudoAbsence",
        Effect == "algorithmGLM:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set2var:algorithmEnsemble" ~ "VarSet*Algo",
        Effect == "algorithmMAXENT.Phillips:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set2var:algorithmMAXENT.Background" ~ "VarSet*Algo",
        Effect == "algorithmMARS:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set2var:algorithmMAXENT.Phillips" ~ "VarSet*Algo",
        Effect == "algorithmMAXENT.Background:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEnsemble:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmMAXENT.Phillips:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmGAM:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmMAXENT.Background:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:algorithmEnsemble" ~ "VarSet*Algo",
        Effect == "algorithmEnsemble:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_setallvar:pseudoabs_setPA=1k" ~ "VarSet*PseudoAbsence",
        #
        Effect == "predictor_setallvar:algorithmEnsemble" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmMAXENT.Background" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmMAXENT.Phillips" ~ "VarSet*Algo",
        Effect == "algorithmGBM:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_setallvar:algorithmMAXENT.Phillips" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmMAXENT.Background" ~ "VarSet*Algo",
        Effect == "algorithmMAXENT.Background:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set2var:pseudoabs_setPA=1k" ~ "VarSet*PseudoAbsence",
        #
        Effect == "predictor_set4var:algorithmEMwmeanByTSS" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmEMmeanByCSI" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmEMmeanByTSS" ~ "VarSet*Algo",
        Effect == "algorithmEMwmeanByTSS" ~ "Algo",
        Effect == "algorithmEMwmeanByCSI:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEMmeanByCSI:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEMmeanByTSS:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:algorithmEMmeanByCSI" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmEMmeanByTSS" ~ "VarSet*Algo",
        Effect == "algorithmEMmeanByTSS:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set4var:algorithmEMwmeanByCSI" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmEMmeanByTSS" ~ "VarSet*Algo",
        Effect == "predictor_set4var:algorithmEMwmeanByCSI" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmEMwmeanByCSI" ~ "VarSet*Algo",
        Effect == "algorithmEMwmeanByTSS:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "predictor_setallvar:algorithmEMwmeanByTSS" ~ "VarSet*Algo",
        Effect == "predictor_set2var:algorithmEMmeanByTSS" ~ "VarSet*Algo",
        Effect == "predictor_set2var:algorithmEMmeanByCSI" ~ "VarSet*Algo",
        Effect == "predictor_set2var:algorithmEMwmeanByCSI" ~ "VarSet*Algo",
        Effect == "algorithmEMwmeanByCSI" ~ "Algo",
        Effect == "algorithmEMwmeanByCSI:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "predictor_set2var:algorithmEMwmeanByTSS" ~ "VarSet*Algo",
        Effect == "algorithmEMmeanByCSI:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEMmeanByTSS:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEMwmeanByTSS:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEMmeanByTSS" ~ "Algo",
        Effect == "algorithmEMmeanByCSI" ~ "Algo",
        Effect == "algorithmMAXENT.Background:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmMAXENT.Phillips:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmEMwmeanByCSI:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmEMmeanByCSI:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmEMmeanByTSS:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmEMwmeanByTSS:niche_overlap_value" ~ "Algo*Niche",
        Effect == "algorithmEM:niche_overlap_value" ~ "Algo*Niche",
        Effect == "predictor_set4var:algorithmEM" ~ "VarSet*Algo",
        Effect == "predictor_set2var:algorithmEM" ~ "VarSet*Algo",
        Effect == "predictor_setallvar:algorithmEM" ~ "VarSet*Algo",
        Effect == "algorithmEM" ~ "Algo",
        #
        Effect == "algorithmEM:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence",
        Effect == "algorithmEM:pseudoabs_setPA=1k" ~ "Algo*PseudoAbsence",
        Effect == "algorithmGAM:pseudoabs_setPA=P" ~ "Algo*PseudoAbsence"),
      r_squared = round((Rsq / sum(Rsq)) * 100, 2),
      rsum = sum(Rsq)
    ) %>% 
    group_by(Effect_grouped) %>%
    mutate(
      r_squared_group = sum(Rsq) / rsum * 100
    ) %>% 
    ungroup()
}
