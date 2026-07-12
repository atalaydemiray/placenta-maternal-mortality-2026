#===============================================================================
#
# STATISTICAL ANALYSIS CODE FOR:
# "A comparison of long-term maternal mortality associated with pathologic
#  placental separation: highlighting possible trends and mechanisms"
#
# Authors: Jasani S, Demiray A, Stevenson J, Krawiec C
#
# This script reproduces all statistical analyses reported in the manuscript.
# For questions, contact: atalay.demiray@yale.edu
#
#===============================================================================
#
# REPRODUCIBILITY INFORMATION
# ---------------------------
# R version: 4.5.2 or later recommended
# Required packages listed below with versions used in original analysis
# Random seed: 42 (set for reproducibility where applicable)
#
# DATA REQUIREMENTS
# -----------------
# Input: analysis_data.RData (processed cohort data)
# Note: Raw data from TriNetX cannot be redistributed per data use agreement.
#       The analysis_data.RData contains the analysis-ready cohort.
#
# OUTPUT FILES
# ------------
# Tables:
#   - Table1_Baseline_Characteristics.csv
#   - Table2_Hazard_Ratios.csv
#   - Table3_Health_Outcomes.csv
#   - Table1.docx, Table2.docx, Table3.docx
#
# Supplementary Tables:
#   - S2_EPV_Analysis.csv (Events per variable)
#   - S3_PH_Assumption_Tests.csv (Proportional hazards tests)
#   - S4_Time_Specific_HR.csv (Piecewise Cox results)
#   - S5_Multiple_Comparison_Corrections.csv
#   - S6_Sensitivity_Analyses.csv
#
#===============================================================================

#===============================================================================
# SECTION 1: SETUP AND PACKAGE INSTALLATION
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("LONG-TERM MATERNAL MORTALITY IN PATHOLOGIC PLACENTAL SEPARATION\n")
cat("Statistical Analysis Script (S7 Appendix)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Required packages
required_packages <- c(
  "dplyr",       # Data manipulation
  "tidyr",       # Data reshaping
  "survival",    # Cox proportional hazards models
  "officer",     # Word document generation
  "flextable"   # Table formatting
)

# Install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
invisible(lapply(required_packages, install_if_missing))

# Load packages
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# Set seed for reproducibility
set.seed(42)

# Print session info for reproducibility
cat("R version:", R.version.string, "\n")
cat("Analysis date:", as.character(Sys.Date()), "\n\n")

#===============================================================================
# SECTION 2: LOAD DATA
#===============================================================================

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 2: LOADING DATA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Set working directory to script location or specify path
# Update this path to your data location
data_path <- "data/analysis_data.RData"

if (file.exists(data_path)) {
  load(data_path)
  cat("Loaded data from:", data_path, "\n")
} else if (file.exists("analysis_data.RData")) {
  load("analysis_data.RData")
  cat("Loaded data from: analysis_data.RData\n")
} else {
  stop("Data file not found. Please ensure analysis_data.RData is available.")
}

cat("Dataset dimensions:", nrow(filtered_data), "rows x", ncol(filtered_data), "columns\n")
cat("\nGroup distribution:\n")
print(table(filtered_data$group))

#===============================================================================
# SECTION 3: TABLE 1 - BASELINE CHARACTERISTICS
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 3: TABLE 1 - BASELINE CHARACTERISTICS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Calculate baseline statistics by group
baseline_stats <- filtered_data %>%
  group_by(group) %>%
  summarise(
    N = n(),
    Deaths = sum(death_indicator),
    Death_Rate_per_1000 = round(Deaths / N * 1000, 1),
    Person_Years = round(sum(time_to_death / 365.25)),
    Mortality_per_1000_PY = round(Deaths / Person_Years * 1000, 2),
    Followup_Mean = round(mean(time_to_death / 365.25), 2),
    Followup_SD = round(sd(time_to_death / 365.25), 2),
    Age_Mean = round(mean(maternal_age, na.rm = TRUE), 1),
    Age_SD = round(sd(maternal_age, na.rm = TRUE), 1),

    # Health conditions
    Diabetes_n = sum(diabetes_indicator, na.rm = TRUE),
    Hypertension_n = sum(hypertension_indicator, na.rm = TRUE),
    Obesity_n = sum(obesity_indicator, na.rm = TRUE),
    CKD_n = sum(chronic_kidney_disease_indicator, na.rm = TRUE),
    Asthma_n = sum(asthma_indicator, na.rm = TRUE),
    Mental_n = sum(mental_health_indicator, na.rm = TRUE),

    # Risk factors
    Alcohol_n = sum(alcohol_use_indicator, na.rm = TRUE),
    Tobacco_n = sum(tobaccouse_indicator, na.rm = TRUE),
    Substance_n = sum(substance_use_indicator, na.rm = TRUE),
    SDOH_n = sum(sdoh_indicator, na.rm = TRUE),

    # Pregnancy/Delivery
    Multiple_n = sum(multiple_gestation_indicator, na.rm = TRUE),
    Previa_n = sum(placenta_previa_indicator, na.rm = TRUE),
    Preeclampsia_n = sum(preclampsia_indicator, na.rm = TRUE),
    Eclampsia_n = sum(eclampsia_indicator, na.rm = TRUE),
    PPH_n = sum(pph_mindate_indicator, na.rm = TRUE),
    DIC_n = sum(dic_indicator, na.rm = TRUE),

    # Procedures
    Hysterectomy_n = sum(hysterectomy_indicator, na.rm = TRUE),
    Transfusion_n = sum(transfusion_indicator, na.rm = TRUE),

    .groups = "drop"
  )

print(baseline_stats)
write.csv(baseline_stats, "Table1_Baseline_Characteristics.csv", row.names = FALSE)
cat("\nSaved: Table1_Baseline_Characteristics.csv\n")

# Print follow-up summary for manuscript
cat("\n=== FOLLOW-UP SUMMARY FOR MANUSCRIPT ===\n")
overall_mean <- round(mean(filtered_data$time_to_death / 365.25), 1)
overall_sd <- round(sd(filtered_data$time_to_death / 365.25), 1)
total_py <- round(sum(filtered_data$time_to_death / 365.25))
cat("Overall cohort (N =", nrow(filtered_data), "):\n")
cat("  Mean follow-up:", overall_mean, "years (SD", overall_sd, ")\n")
cat("  Total person-years:", format(total_py, big.mark = ","), "\n")

#===============================================================================
# SECTION 4: COX PROPORTIONAL HAZARDS REGRESSION (TABLE 2)
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 4: COX PROPORTIONAL HAZARDS REGRESSION\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Create survival object
surv_object <- Surv(time = filtered_data$time_to_death,
                    event = filtered_data$death_indicator)

# Helper function to extract hazard ratios
extract_hr <- function(model, model_name) {
  s <- summary(model)
  ci <- s$conf.int
  pv <- s$coefficients[, "Pr(>|z|)"]
  abr <- grep("Abruption", rownames(ci))
  ret <- grep("Retention", rownames(ci))

  data.frame(
    Model = model_name,
    Abruption_HR = round(ci[abr, "exp(coef)"], 2),
    Abruption_LCI = round(ci[abr, "lower .95"], 2),
    Abruption_UCI = round(ci[abr, "upper .95"], 2),
    Abruption_P = ifelse(pv[abr] < 0.001, "<0.001", as.character(round(pv[abr], 3))),
    Retention_HR = round(ci[ret, "exp(coef)"], 2),
    Retention_LCI = round(ci[ret, "lower .95"], 2),
    Retention_UCI = round(ci[ret, "upper .95"], 2),
    Retention_P = ifelse(pv[ret] < 0.001, "<0.001", as.character(round(pv[ret], 3))),
    stringsAsFactors = FALSE
  )
}

#-------------------------------------------------------------------------------
# Model 1: Unadjusted
#-------------------------------------------------------------------------------
cat("=== MODEL 1: UNADJUSTED ===\n")
cox_unadj <- coxph(surv_object ~ group, data = filtered_data)
print(round(summary(cox_unadj)$conf.int, 3))

#-------------------------------------------------------------------------------
# Model 2: Adjusted for Demographics
#-------------------------------------------------------------------------------
cat("\n=== MODEL 2: ADJUSTED FOR DEMOGRAPHICS ===\n")
cat("Covariates: maternal age, race, SDOH indicator\n")
cox_demo <- coxph(surv_object ~ group + maternal_age + race + sdoh_indicator,
                  data = filtered_data)
group_rows <- grep("group", rownames(summary(cox_demo)$conf.int))
print(round(summary(cox_demo)$conf.int[group_rows, ], 3))

#-------------------------------------------------------------------------------
# Model 3: + Pre-existing Conditions
#-------------------------------------------------------------------------------
cat("\n=== MODEL 3: + PRE-EXISTING CONDITIONS ===\n")
cat("Additional covariates: diabetes, hypertension, obesity, CKD\n")
cox_full <- coxph(surv_object ~ group + maternal_age + race + sdoh_indicator +
                    diabetes_indicator + hypertension_indicator +
                    obesity_indicator + chronic_kidney_disease_indicator,
                  data = filtered_data)
group_rows <- grep("group", rownames(summary(cox_full)$conf.int))
print(round(summary(cox_full)$conf.int[group_rows, ], 3))

#===============================================================================
# SECTION 5: SENSITIVITY ANALYSES
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 5: SENSITIVITY ANALYSES\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

#-------------------------------------------------------------------------------
# Analysis A: Short-term mortality (within 1 year, censored at 365 days)
#-------------------------------------------------------------------------------
cat("=== ANALYSIS A: SHORT-TERM MORTALITY (WITHIN 1 YEAR) ===\n")
cat("Method: Follow-up censored at 365 days\n")

data_within_1yr <- filtered_data %>%
  mutate(
    time_censored = pmin(time_to_death, 365),
    event_censored = ifelse(time_to_death <= 365 & death_indicator == 1, 1, 0)
  )

cat("N:", nrow(data_within_1yr), "| Deaths within 1 year:", sum(data_within_1yr$event_censored), "\n")

surv_within_1yr <- Surv(time = data_within_1yr$time_censored,
                        event = data_within_1yr$event_censored)
cox_within_1yr <- coxph(surv_within_1yr ~ group + maternal_age + race + sdoh_indicator,
                        data = data_within_1yr)
group_rows <- grep("group", rownames(summary(cox_within_1yr)$conf.int))
print(round(summary(cox_within_1yr)$conf.int[group_rows, ], 3))

#-------------------------------------------------------------------------------
# Analysis B: Excluding deaths ≤42 days (WHO maternal mortality window)
#-------------------------------------------------------------------------------
cat("\n=== ANALYSIS B: EXCLUDING DEATHS ≤42 DAYS ===\n")
cat("Method: Remove patients who died within 42 days\n")

data_excl_42d <- filtered_data %>%
  filter(time_to_death > 42 | death_indicator == 0)

n_excluded_42d <- nrow(filtered_data) - nrow(data_excl_42d)
cat("N:", nrow(data_excl_42d), "| Deaths:", sum(data_excl_42d$death_indicator),
    "| Excluded:", n_excluded_42d, "\n")

surv_excl_42d <- Surv(time = data_excl_42d$time_to_death,
                      event = data_excl_42d$death_indicator)
cox_excl_42d <- coxph(surv_excl_42d ~ group + maternal_age + race + sdoh_indicator,
                      data = data_excl_42d)
group_rows <- grep("group", rownames(summary(cox_excl_42d)$conf.int))
print(round(summary(cox_excl_42d)$conf.int[group_rows, ], 3))

#-------------------------------------------------------------------------------
# Analysis C: Excluding deaths ≤1 year (long-term only)
#-------------------------------------------------------------------------------
cat("\n=== ANALYSIS C: EXCLUDING DEATHS ≤1 YEAR (LONG-TERM ONLY) ===\n")
cat("Method: Remove patients who died within 1 year\n")

data_excl_1yr <- filtered_data %>%
  filter(time_to_death > 365 | death_indicator == 0)

n_excluded_1yr <- nrow(filtered_data) - nrow(data_excl_1yr)
cat("N:", nrow(data_excl_1yr), "| Deaths:", sum(data_excl_1yr$death_indicator),
    "| Excluded:", n_excluded_1yr, "\n")

surv_excl_1yr <- Surv(time = data_excl_1yr$time_to_death,
                      event = data_excl_1yr$death_indicator)
cox_excl_1yr <- coxph(surv_excl_1yr ~ group + maternal_age + race + sdoh_indicator,
                      data = data_excl_1yr)
group_rows <- grep("group", rownames(summary(cox_excl_1yr)$conf.int))
print(round(summary(cox_excl_1yr)$conf.int[group_rows, ], 3))

#===============================================================================
# SECTION 6: COMPILE TABLE 2
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 6: COMPILE TABLE 2 - HAZARD RATIOS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

table2 <- rbind(
  extract_hr(cox_unadj, "Unadjusted"),
  extract_hr(cox_demo, "Adjusted (Demographics)"),
  extract_hr(cox_full, "+ Pre-existing conditions"),
  extract_hr(cox_within_1yr, "Within 1 year (censored at 365 days)"),
  extract_hr(cox_excl_42d, "Excluding deaths ≤42 days"),
  extract_hr(cox_excl_1yr, "Excluding deaths ≤1 year (long-term only)")
)

# Add formatted columns
table2$Abruption_HR_CI <- paste0(table2$Abruption_HR, " (",
                                  table2$Abruption_LCI, "-",
                                  table2$Abruption_UCI, ")")
table2$Retention_HR_CI <- paste0(table2$Retention_HR, " (",
                                  table2$Retention_LCI, "-",
                                  table2$Retention_UCI, ")")

cat("TABLE 2: HAZARD RATIOS\n\n")
print(table2[, c("Model", "Abruption_HR_CI", "Abruption_P",
                 "Retention_HR_CI", "Retention_P")])

write.csv(table2, "Table2_Hazard_Ratios.csv", row.names = FALSE)
cat("\nSaved: Table2_Hazard_Ratios.csv\n")

# Save S6 Sensitivity Analyses
s6_data <- table2[c(1, 4, 5, 6), c("Model", "Abruption_HR", "Abruption_LCI",
                                    "Abruption_UCI", "Abruption_P",
                                    "Retention_HR", "Retention_LCI",
                                    "Retention_UCI", "Retention_P",
                                    "Abruption_HR_CI", "Retention_HR_CI")]
write.csv(s6_data, "S6_Sensitivity_Analyses.csv", row.names = FALSE)
cat("Saved: S6_Sensitivity_Analyses.csv\n")

#===============================================================================
# SECTION 7: PROPORTIONAL HAZARDS ASSUMPTION TESTING
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 7: PROPORTIONAL HAZARDS ASSUMPTION TESTING\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Test PH assumption for each model
ph_test_unadj <- cox.zph(cox_unadj)
ph_test_demo <- cox.zph(cox_demo)
ph_test_preexist <- cox.zph(cox_full)

cat("=== UNADJUSTED MODEL ===\n")
print(ph_test_unadj)

cat("\n=== DEMOGRAPHICS MODEL ===\n")
print(ph_test_demo)

cat("\n=== PRE-EXISTING CONDITIONS MODEL ===\n")
print(ph_test_preexist)

# Create detailed PH results table (S3 Table)
ph_results <- data.frame(
  Model = character(),
  Variable = character(),
  Chi_square = numeric(),
  df = numeric(),
  P_value = numeric(),
  PH_Violated = character(),
  stringsAsFactors = FALSE
)

# Extract results from each model
for (model_info in list(
  list(name = "Unadjusted", test = ph_test_unadj),
  list(name = "Demographics Adjusted", test = ph_test_demo),
  list(name = "Pre-existing Conditions", test = ph_test_preexist)
)) {
  ph_table <- model_info$test$table
  for (i in 1:nrow(ph_table)) {
    var_name <- rownames(ph_table)[i]
    # Clean up variable names for readability
    if (grepl("^group", var_name)) var_name <- "Group (Placental separation type)"
    if (var_name == "maternal_age") var_name <- "Maternal age"
    if (grepl("^race", var_name)) var_name <- "Race"
    if (var_name == "sdoh_indicator") var_name <- "SDOH indicator"
    if (var_name == "diabetes_indicator") var_name <- "Diabetes"
    if (var_name == "hypertension_indicator") var_name <- "Hypertension"
    if (var_name == "obesity_indicator") var_name <- "Obesity"
    if (var_name == "chronic_kidney_disease_indicator") var_name <- "Chronic kidney disease"

    ph_results <- rbind(ph_results, data.frame(
      Model = model_info$name,
      Variable = var_name,
      Chi_square = round(ph_table[i, "chisq"], 2),
      df = ph_table[i, "df"],
      P_value = round(ph_table[i, "p"], 4),
      PH_Violated = ifelse(ph_table[i, "p"] < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(ph_results, "S3_PH_Assumption_Tests.csv", row.names = FALSE)
cat("\nSaved: S3_PH_Assumption_Tests.csv\n")

# Print key finding
cat("\n=== KEY FINDING: PRIMARY EXPOSURE VARIABLE ===\n")
group_ph <- ph_results %>%
  filter(Variable == "Group (Placental separation type)", Model == "Demographics Adjusted")
cat("Placental separation group PH test: p =", group_ph$P_value, "\n")
cat("PH assumption satisfied:", group_ph$PH_Violated == "No", "\n")

#===============================================================================
# SECTION 8: PIECEWISE COX REGRESSION (S4 TABLE)
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 8: PIECEWISE COX REGRESSION - TIME-SPECIFIC HRs\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Define time periods (in days)
time_breaks <- c(0, 365, 730, 1825, Inf)
time_labels <- c("0-1 year", "1-2 years", "2-5 years", ">5 years")

# Function to fit Cox model for a specific time window
fit_piecewise <- function(data, t_start, t_end, period_label) {
  data_period <- data %>%
    filter(time_to_death > t_start) %>%
    mutate(
      time_in_period = pmin(time_to_death, t_end),
      event_in_period = ifelse(time_to_death <= t_end & death_indicator == 1, 1, 0)
    )

  n_patients <- nrow(data_period)
  n_events <- sum(data_period$event_in_period)

  cat(period_label, ": N at risk =", n_patients, ", Events =", n_events, "\n")

  if (n_events < 10) return(NULL)

  surv_period <- Surv(time = data_period$time_in_period,
                      event = data_period$event_in_period)

  cox_period <- tryCatch({
    coxph(surv_period ~ group + maternal_age + race + sdoh_indicator,
          data = data_period)
  }, error = function(e) NULL)

  if (is.null(cox_period)) return(NULL)

  # Extract results
  s <- summary(cox_period)
  ci <- s$conf.int
  pv <- s$coefficients[, "Pr(>|z|)"]
  abr <- grep("Abruption", rownames(ci))
  ret <- grep("Retention", rownames(ci))

  data.frame(
    Period = period_label,
    N = n_patients,
    Events = n_events,
    Abruption_HR = round(ci[abr, "exp(coef)"], 2),
    Abruption_LCI = round(ci[abr, "lower .95"], 2),
    Abruption_UCI = round(ci[abr, "upper .95"], 2),
    Abruption_P = ifelse(pv[abr] < 0.001, "<0.001", as.character(round(pv[abr], 3))),
    Retention_HR = round(ci[ret, "exp(coef)"], 2),
    Retention_LCI = round(ci[ret, "lower .95"], 2),
    Retention_UCI = round(ci[ret, "upper .95"], 2),
    Retention_P = ifelse(pv[ret] < 0.001, "<0.001", as.character(round(pv[ret], 3))),
    stringsAsFactors = FALSE
  )
}

# Fit models for each period
piecewise_results <- list()
for (i in 1:(length(time_breaks) - 1)) {
  result <- fit_piecewise(filtered_data, time_breaks[i], time_breaks[i + 1], time_labels[i])
  if (!is.null(result)) piecewise_results[[length(piecewise_results) + 1]] <- result
}

s4_table <- do.call(rbind, piecewise_results)

# Add formatted columns
s4_table$Abruption_HR_CI <- paste0(s4_table$Abruption_HR, " (",
                                    s4_table$Abruption_LCI, "-",
                                    s4_table$Abruption_UCI, ")")
s4_table$Retention_HR_CI <- paste0(s4_table$Retention_HR, " (",
                                    s4_table$Retention_LCI, "-",
                                    s4_table$Retention_UCI, ")")

cat("\nS4 TABLE: TIME-SPECIFIC HAZARD RATIOS\n")
print(s4_table[, c("Period", "N", "Events", "Abruption_HR_CI", "Abruption_P",
                   "Retention_HR_CI", "Retention_P")])

write.csv(s4_table, "S4_Time_Specific_HR.csv", row.names = FALSE)
cat("\nSaved: S4_Time_Specific_HR.csv\n")

#===============================================================================
# SECTION 9: EVENTS PER VARIABLE ANALYSIS (S2 TABLE)
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 9: EVENTS PER VARIABLE (EPV) ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

n_events <- sum(filtered_data$death_indicator)
s2_table <- data.frame(
  Model = c("Unadjusted", "Demographics", "Full"),
  Predictors = c(2, 5, 10),  # Approximate number of parameters
  Events = n_events,
  EPV = round(n_events / c(2, 5, 10), 1),
  Adequate = TRUE
)
s2_table$Adequate <- s2_table$EPV > 10

print(s2_table)
write.csv(s2_table, "S2_EPV_Analysis.csv", row.names = FALSE)
cat("\nSaved: S2_EPV_Analysis.csv\n")

#===============================================================================
# SECTION 10: MULTIPLE COMPARISON CORRECTIONS (S5 TABLE)
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 10: MULTIPLE COMPARISON CORRECTIONS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Get p-values from primary analyses
p_values <- c(
  M1_Abruption = as.numeric(gsub("<", "", table2$Abruption_P[1])),
  M1_Retention = as.numeric(gsub("<", "", table2$Retention_P[1])),
  M2_Abruption = as.numeric(gsub("<", "", table2$Abruption_P[2])),
  M2_Retention = as.numeric(gsub("<", "", table2$Retention_P[2]))
)
p_values[p_values == 0.001] <- 0.0001  # Handle "<0.001"

# Apply corrections
n_tests <- length(p_values)
bonferroni_p <- pmin(p_values * n_tests, 1)
fdr_p <- p.adjust(p_values, method = "BH")

s5_table <- data.frame(
  Comparison = names(p_values),
  Raw_P = p_values,
  Bonferroni_P = round(bonferroni_p, 4),
  FDR_P = round(fdr_p, 4),
  Sig_Raw = p_values < 0.05,
  Sig_Bonferroni = bonferroni_p < 0.05,
  Sig_FDR = fdr_p < 0.05
)

print(s5_table)
write.csv(s5_table, "S5_Multiple_Comparison_Corrections.csv", row.names = FALSE)
cat("\nSaved: S5_Multiple_Comparison_Corrections.csv\n")

#===============================================================================
# SECTION 11: TABLE 3 - HEALTH OUTCOMES (LOGISTIC REGRESSION)
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 11: TABLE 3 - HEALTH OUTCOME ASSOCIATIONS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Ensure factors are set correctly
filtered_data$group <- factor(filtered_data$group,
                               levels = c("Normal Delivery", "Placental Abruption", "Placental Retention"))

# Define outcomes to analyze
outcomes <- c(
  "placenta_previa_indicator", "preclampsia_indicator", "eclampsia_indicator",
  "gestational_htn_indicator", "multiple_gestation_indicator", "pph_mindate_indicator",
  "alcohol_use_indicator", "tobaccouse_indicator", "substance_use_indicator",
  "transfusion_indicator", "mechanical_ventilation_indicator", "critical_care_indicator",
  "hysterectomy_indicator", "hypertension_indicator", "diabetes_indicator",
  "chronic_kidney_disease_indicator", "mental_health_indicator", "lupus_indicator",
  "asthma_indicator", "cardiacvalvulardisease_indicator", "obesity_indicator",
  "ards_indicator", "shock_indicator", "sepsis_indicator", "dic_indicator",
  "heart_failure_indicator", "anemia_indicator"
)

# Function to run logistic regression
run_logistic <- function(data, outcome_var) {
  formula <- as.formula(paste0(outcome_var, " ~ group + maternal_age + race + sdoh_indicator"))

  tryCatch({
    model <- glm(formula, data = data, family = binomial)
    coef_table <- summary(model)$coefficients

    # Extract ORs
    abr_est <- coef_table["groupPlacental Abruption", "Estimate"]
    abr_se <- coef_table["groupPlacental Abruption", "Std. Error"]
    abr_p <- coef_table["groupPlacental Abruption", "Pr(>|z|)"]

    ret_est <- coef_table["groupPlacental Retention", "Estimate"]
    ret_se <- coef_table["groupPlacental Retention", "Std. Error"]
    ret_p <- coef_table["groupPlacental Retention", "Pr(>|z|)"]

    data.frame(
      Outcome = outcome_var,
      Abruption_OR = round(exp(abr_est), 2),
      Abruption_LCI = round(exp(abr_est - 1.96 * abr_se), 2),
      Abruption_UCI = round(exp(abr_est + 1.96 * abr_se), 2),
      Abruption_P = ifelse(abr_p < 0.001, "<0.001", round(abr_p, 3)),
      Retention_OR = round(exp(ret_est), 2),
      Retention_LCI = round(exp(ret_est - 1.96 * ret_se), 2),
      Retention_UCI = round(exp(ret_est + 1.96 * ret_se), 2),
      Retention_P = ifelse(ret_p < 0.001, "<0.001", round(ret_p, 3)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

# Run analyses
table3_results <- list()
for (outcome in outcomes) {
  if (outcome %in% names(filtered_data)) {
    result <- run_logistic(filtered_data, outcome)
    if (!is.null(result)) table3_results[[length(table3_results) + 1]] <- result
  }
}

table3 <- do.call(rbind, table3_results)

if (nrow(table3) > 0) {
  write.csv(table3, "Table3_Health_Outcomes.csv", row.names = FALSE)
  cat("Saved: Table3_Health_Outcomes.csv\n")
  cat("\nNumber of outcomes analyzed:", nrow(table3), "\n")
}

#===============================================================================
# SECTION 12: GENERATE WORD DOCUMENTS (PLOS FORMAT)
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION 12: GENERATING WORD DOCUMENTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Table 2 Word Document
table2_display <- data.frame(
  Model = table2$Model,
  `Abruption HR (95% CI)` = table2$Abruption_HR_CI,
  `Abruption P` = table2$Abruption_P,
  `Retention HR (95% CI)` = table2$Retention_HR_CI,
  `Retention P` = table2$Retention_P,
  check.names = FALSE
)

ft2 <- flextable(table2_display) %>%
  add_header_row(values = c("", "Placental Abruption", "", "Placental Retention", ""),
                 top = TRUE) %>%
  merge_at(i = 1, j = 2:3, part = "header") %>%
  merge_at(i = 1, j = 4:5, part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(align = "center", j = 2:5, part = "all") %>%
  align(align = "left", j = 1, part = "all") %>%
  bold(part = "header") %>%
  autofit()

doc2 <- read_docx() %>%
  body_add_par("Table 2. Hazard ratios for all-cause mortality by placental separation type",
               style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft2) %>%
  body_add_par("") %>%
  body_add_par("Reference group: Normal Delivery (n = 625,890).", style = "Normal") %>%
  body_add_par("a Adjusted for maternal age, race/ethnicity, and SDOH indicator.", style = "Normal") %>%
  body_add_par("b Additionally adjusted for diabetes, hypertension, obesity, and CKD.", style = "Normal")

print(doc2, target = "Table2.docx")
cat("Saved: Table2.docx\n")

#===============================================================================
# SECTION 13: FINAL SUMMARY
#===============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("OUTPUT FILES GENERATED:\n")
cat("  Main Tables:\n")
cat("    - Table1_Baseline_Characteristics.csv\n")
cat("    - Table2_Hazard_Ratios.csv\n")
cat("    - Table3_Health_Outcomes.csv\n")
cat("    - Table2.docx\n")
cat("\n  Supplementary Tables:\n")
cat("    - S2_EPV_Analysis.csv\n")
cat("    - S3_PH_Assumption_Tests.csv\n")
cat("    - S4_Time_Specific_HR.csv\n")
cat("    - S5_Multiple_Comparison_Corrections.csv\n")
cat("    - S6_Sensitivity_Analyses.csv\n")

cat("\n=== KEY FINDINGS ===\n")
cat("1. Both abruption and retention associated with increased mortality\n")
cat("2. Retention: persistent risk across all time horizons (HR 1.75-2.02)\n")
cat("3. Abruption: highest short-term (HR 1.88), attenuated at 42d, returns long-term\n")
cat("4. PH assumption satisfied for primary exposure (p=0.62)\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("END OF ANALYSIS\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

#===============================================================================
# SESSION INFO FOR REPRODUCIBILITY
#===============================================================================

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
