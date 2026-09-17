# 1. LOAD REQUIRED PACKAGES
library(dplyr)     # Data wrangling
library(car)       # Anova and VIF testing
library(agricolae) # Post-hoc LSD/Tukey tests
library(corrplot)  # Correlation matrix visualization
library(MASS)      # Stepwise regression (stepAIC)
library(vegan)     # Redundancy Analysis (RDA)
library(ggplot2)   # Publication graphics
library(ggpubr)    # Multi-panel plots
library(Hmisc)     # Pearson correlation matrix with p-values
library(tibble)    # Data frame utilities

# 2. DATA IMPORT & CLEANING
df <- read.csv("C:/Users/user/OneDrive/Desktop/thesis complete data.csv")  
colnames(df) <- make.names(colnames(df))  

# Clean and standardize variable definitions globally
df <- df %>%  
  mutate( 
    Land_Type = as.factor(Land.Type),  
    pH        = as.numeric(pH),  
    Clay      = as.numeric(Clay....),  
    OM        = as.numeric(OM....),  
    N         = as.numeric(N....),  
    P2O5      = as.numeric(P2O5.kg.ha.),  
    K         = as.numeric(K..kg.ha.),  
    S         = as.numeric(S),  
    k         = as.numeric(k..day.1.),
    Mr_g      = as.numeric(Mr.g),
    Mr_r      = as.numeric(Mr.r)
  ) 

# 3. SUMMARY STATISTICS BY LAND USE
summary_stats <- df %>%  
  group_by(Land_Type) %>%  
  summarise( 
    pH_mean = mean(pH, na.rm=TRUE),   pH_sd = sd(pH, na.rm=TRUE),  
    OM_mean = mean(OM, na.rm=TRUE),   OM_sd = sd(OM, na.rm=TRUE),  
    N_mean  = mean(N, na.rm=TRUE),    N_sd  = sd(N, na.rm=TRUE),  
    P_mean  = mean(P2O5, na.rm=TRUE), P_sd  = sd(P2O5, na.rm=TRUE),  
    K_mean  = mean(K, na.rm=TRUE),    K_sd  = sd(K, na.rm=TRUE),  
    S_mean  = mean(S, na.rm=TRUE),    S_sd  = sd(S, na.rm=TRUE),  
    k_mean  = mean(k, na.rm=TRUE),    k_sd  = sd(k, na.rm=TRUE),
    .groups = "drop"
  )  
cat("=== SUMMARY STATISTICS ===\n")  
print(summary_stats) 
# SECTION 3.5: DIAGNOSTIC EVALUATION OF DATA NORMALITY (SHAPIRO-WILK TESTS)
# ==============================================================================
cat("\n=================================================================\n")
cat(" SECTION 3.5: SHAPIRO-WILK NORMALITY TESTS \n")
cat("=================================================================\n\n")

all_vars <- c("pH", "Clay", "OM", "N", "P2O5", "K", "Mr_g", "Mr_r", "S", "k")

# 1. Overall Dataset Normality (N = 30)
cat("--- 1. OVERALL DATASET NORMALITY (N = 30) ---\n")
overall_norm <- data.frame(Variable = character(), W = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (v in all_vars) {
  sw <- shapiro.test(df[[v]])
  overall_norm <- rbind(overall_norm, data.frame(Variable = v, W = round(sw$statistic, 4), p_value = round(sw$p.value, 4)))
}
print(overall_norm, row.names = FALSE)

# 2. Within-Group Normality Range by Land Type (n = 10 each)
cat("\n--- 2. WITHIN-GROUP NORMALITY RANGE BY LAND TYPE (n = 10) ---\n")
group_norm <- data.frame(Variable = character(), Min_p = numeric(), Max_p = numeric(), Range_Fmt = character(), stringsAsFactors = FALSE)

for (v in all_vars) {
  p_vals <- df %>%
    group_by(Land_Type) %>%
    summarise(p = shapiro.test(.data[[v]])$p.value, .groups = "drop") %>%
    pull(p)
  
  min_p <- min(p_vals)
  max_p <- max(p_vals)
  
  group_norm <- rbind(group_norm, data.frame(
    Variable = v,
    Min_p = round(min_p, 4),
    Max_p = round(max_p, 4),
    Range_Fmt = sprintf("%.4f - %.4f", min_p, max_p)
  ))
}
print(group_norm[, c("Variable", "Range_Fmt")], row.names = FALSE)

# 3. OLS Regression Residual Normality
cat("\n--- 3. OLS REGRESSION RESIDUAL NORMALITY ---\n")
full_k_model <- lm(k ~ OM + N + P2O5 + K + pH + Clay, data = df)
full_S_model <- lm(S ~ OM + N + P2O5 + K + pH + Clay, data = df)

sw_k_res <- shapiro.test(residuals(full_k_model))
sw_S_res <- shapiro.test(residuals(full_S_model))

cat(sprintf("Decomposition Rate (k) Residuals: W = %.4f | p-value = %.4f\n", sw_k_res$statistic, sw_k_res$p.value))
cat(sprintf("Stabilization Factor (S) Residuals: W = %.4f | p-value = %.4f\n", sw_S_res$statistic, sw_S_res$p.value))
cat("=================================================================\n\n")

# SECTION 4.1: ANOVA & TUKEY HSD FOR SOIL PROPERTIES
vars_sec_4_1 <- c("pH", "Clay", "OM", "N", "P2O5", "K") 
cat("\n=================================================================\n")  
cat(" EXACT METRICS FOR CHAPTER 4.1 & TABLE 4.1 \n")  
cat("=================================================================\n\n")

for (v in vars_sec_4_1) { 
  fit <- aov(as.formula(paste(v, "~ Land_Type")), data = df) 
  tukey_res <- agricolae::HSD.test(fit, "Land_Type", group = TRUE) 
  
  groups <- tukey_res$groups %>%  
    tibble::rownames_to_column(var = "Land_Type") %>%  
    dplyr::rename(Letter = groups) 
  
  stats <- df %>%  
    group_by(Land_Type) %>%  
    summarise( 
      Mean = mean(.data[[v]], na.rm = TRUE),  
      SD   = sd(.data[[v]], na.rm = TRUE),  
      .groups = "drop" 
    ) %>%  
    left_join(groups, by = "Land_Type") %>%  
    mutate(Formatted = sprintf("%.3f ± %.3f (%s)", Mean, SD, Letter)) 
  
  aov_summary <- summary(fit)[[1]]  
  f_val       <- aov_summary$`F value`[1]  
  p_val       <- aov_summary$`Pr(>F)`[1] 
  
  cat(sprintf("=== VARIABLE: %s ===\n", v))  
  cat(sprintf("ANOVA: F(2, 27) = %.3f | p-value = %.4f\n", f_val, p_val))  
  print(stats %>% dplyr::select(Land_Type, Formatted))  
  cat("\n-----------------------------------------------------------------\n")
} 

# SECTION 4.2: ANOVA & TUKEY HSD FOR TBI PARAMETERS
vars_sec_4_2 <- c("Mr_g", "Mr_r", "S", "k") 
cat("\n=================================================================\n")  
cat(" EXACT METRICS FOR CHAPTER 4.2 & TABLE 4.2 \n")  
cat("=================================================================\n\n")

for (v in vars_sec_4_2) { 
  fit <- aov(as.formula(paste(v, "~ Land_Type")), data = df) 
  tukey_res <- agricolae::HSD.test(fit, "Land_Type", group = TRUE) 
  
  groups <- tukey_res$groups %>%  
    tibble::rownames_to_column(var = "Land_Type") %>%  
    dplyr::rename(Letter = groups) 
  
  stats <- df %>%  
    group_by(Land_Type) %>%  
    summarise( 
      Mean = mean(.data[[v]], na.rm = TRUE),  
      SD   = sd(.data[[v]], na.rm = TRUE),  
      SE   = sd(.data[[v]], na.rm = TRUE) / sqrt(n()),  
      .groups = "drop" 
    ) %>%  
    left_join(groups, by = "Land_Type") %>%  
    mutate( 
      Formatted_SD = sprintf("%.5f ± %.5f (%s)", Mean, SD, Letter),  
      Formatted_SE = sprintf("%.5f ± %.5f (%s)", Mean, SE, Letter) 
    ) 
  
  aov_summary <- summary(fit)[[1]]  
  f_val       <- aov_summary$`F value`[1]  
  p_val       <- aov_summary$`Pr(>F)`[1] 
  
  cat(sprintf("=== VARIABLE: %s ===\n", v))  
  cat(sprintf("ANOVA: F(2, 27) = %.3f | p-value = %.6f\n", f_val, p_val))  
  cat("Mean ± SD:\n")  
  print(stats %>% dplyr::select(Land_Type, Formatted_SD))  
  cat("\nMean ± SE:\n")  
  print(stats %>% dplyr::select(Land_Type, Formatted_SE))  
  cat("\n-----------------------------------------------------------------\n")
} 

# SECTION 4.3.1: PEARSON CORRELATION ANALYSIS FOR k AND S
soil_vars <- c("OM", "N", "K", "pH", "P2O5", "Clay") 
cor_data  <- df[, c("k", "S", soil_vars)] 
cor_results <- Hmisc::rcorr(as.matrix(cor_data), type = "pearson") 

pearson_summary <- data.frame( 
  Variable = c("Soil Organic Matter (OM)", "Total Nitrogen (N)",  
               "Extractable Potassium (K)", "Soil Reaction (pH)",  
               "Available P2O5", "Clay Content (%)"), 
  R_Column = soil_vars, 
  r_k      = cor_results$r["k", soil_vars], 
  p_k      = cor_results$P["k", soil_vars], 
  r_S      = cor_results$r["S", soil_vars], 
  p_S      = cor_results$P["S", soil_vars]
) 

pearson_summary <- pearson_summary %>%  
  mutate( 
    p_k_fmt = case_when( 
      p_k < 0.0001 ~ "< 0.0001***", 
      p_k < 0.001  ~ sprintf("%.4f***", p_k), 
      p_k < 0.01   ~ sprintf("%.4f**", p_k), 
      p_k < 0.05   ~ sprintf("%.4f*", p_k), 
      TRUE         ~ sprintf("%.4f (ns)", p_k) 
    ), 
    p_S_fmt = case_when( 
      p_S < 0.0001 ~ "< 0.0001***", 
      p_S < 0.001  ~ sprintf("%.4f***", p_S), 
      p_S < 0.01   ~ sprintf("%.4f**", p_S), 
      p_S < 0.05   ~ sprintf("%.4f*", p_S), 
      TRUE         ~ sprintf("%.4f (ns)", p_S) 
    ) 
  ) 

cat("\n=================================================================\n")  
cat(" PEARSON CORRELATION MATRIX RESULTS (SECTION 4.3.1) \n")  
cat("=================================================================\n\n")
print(pearson_summary[, c("Variable", "r_k", "p_k_fmt", "r_S", "p_S_fmt")], row.names = FALSE) 

# SECTION 4.3.2: OLS MULTIPLE LINEAR REGRESSION MODELS (k & S)
print_ols_metrics <- function(model, model_name) {
  sum_mod <- summary(model)
  f_stat  <- sum_mod$fstatistic[1]
  df1     <- sum_mod$fstatistic[2]
  df2     <- sum_mod$fstatistic[3]
  p_val   <- pf(f_stat, df1, df2, lower.tail = FALSE)
  
  cat("=================================================================\n")  
  cat(sprintf(" %s RESULTS \n", model_name))  
  cat("=================================================================\n\n")  
  cat(sprintf("R-squared:         %.4f\n", sum_mod$r.squared))  
  cat(sprintf("Adjusted R-squared: %.4f\n", sum_mod$adj.r.squared))  
  cat(sprintf("F-statistic:        %.2f on %d and %d DF\n", f_stat, df1, df2))  
  cat(sprintf("Model p-value:     %.6e\n\n", p_val))  
  cat("Coefficients Table:\n")  
  print(sum_mod$coefficients) 
}

# 1. Stepwise OLS Regression Model for Initial Decomposition Rate (k)  
full_k <- lm(k ~ OM + N + P2O5 + K + pH + Clay, data = df) 
step_k <- MASS::stepAIC(full_k, direction = "both", trace = FALSE) 
print_ols_metrics(step_k, "MODEL 1: INITIAL DECOMPOSITION RATE (k)")

# 2. Stepwise OLS Regression Model for Soil Stabilization Factor (S)  
full_S <- lm(S ~ OM + N + P2O5 + K + pH + Clay, data = df) 
step_S <- MASS::stepAIC(full_S, direction = "both", trace = FALSE) 
print_ols_metrics(step_S, "MODEL 2: SOIL STABILIZATION FACTOR (S)")

# SECTION 4.3.3: MULTIVARIATE REDUNDANCY ANALYSIS (RDA) VERIFICATION
Y <- df[, c("k", "S")]  
X <- df[, c("OM", "N", "K", "pH", "P2O5", "Clay")] 

# Standardize response variables for RDA  
Y_std <- decostand(Y, method = "standardize") 

# Fit RDA model  
rda_model <- rda(Y_std ~ OM + N + K + pH + P2O5 + Clay, data = X) 

# Extract inertia partitioning  
tot_inertia   <- rda_model$tot.chi  
con_inertia   <- rda_model$CCA$tot.chi  
uncon_inertia <- rda_model$CA$tot.chi  

prop_constrained      <- (con_inertia / tot_inertia) * 100 
rda1_eigenvalue       <- rda_model$CCA$eig[1] 
rda1_prop_constrained <- (rda1_eigenvalue / con_inertia) * 100 
rda1_prop_total       <- (rda1_eigenvalue / tot_inertia) * 100 

cat("\n=================================================================\n") 
cat(" SECTION 4.3.3 RDA INERTIA & VARIANCE BREAKDOWN \n") 
cat("=================================================================\n\n") 
cat(sprintf("Total Inertia (Variance): %.6f\n", tot_inertia)) 
cat(sprintf("Constrained Inertia:     %.6f (%.2f%% of Total)\n", con_inertia, prop_constrained)) 
cat(sprintf("Unconstrained (Residual):%.6f (%.2f%% of Total)\n\n", uncon_inertia, 100 - prop_constrained)) 

cat("Axis Breakdown:\n") 
cat(sprintf("RDA1 Eigenvalue (Inertia):    %.8f\n", rda1_eigenvalue)) 
cat(sprintf("RDA1 %% of Constrained Inertia: %.2f%%\n", rda1_prop_constrained)) 
cat(sprintf("RDA1 %% of Total Variance:     %.2f%%\n\n", rda1_prop_total)) 

cat("=================================================================\n") 
cat(" PERMUTATION TEST RESULTS (999 PERMUTATIONS) \n") 
cat("=================================================================\n\n")

set.seed(123) 
cat("1. Axis Significance Test:\n") 
print(anova.cca(rda_model, by = "axis", permutations = 999)) 

cat("\n2. Marginal Predictor Term Significance:\n") 
print(anova.cca(rda_model, by = "term", permutations = 999))
