dir.create("data")
dir.create("scripts")
dir.create("results")
dir.create("Figures")
library(survival)
library(survminer)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

# Create query for limited breast cancer data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  sample.type = "Primary Tumor",  # Only tumor samples
                  access = "open",  # Only open access data
                  experimental.strategy = "RNA-Seq")  # Only RNA-seq data

# Check before download
print(query)
GDCdownload(query)
brca_data <- GDCprepare(query)

# Save the prepared data 
saveRDS(brca_data, "data/brca_data.rds")
print("Data successfully saved!")

brca_data <- readRDS('data/brca_data.rds')
print("Data successfully loaded")

dim(brca_data)

clinical_data <- colData(brca_data)
print("Clinical variables available:")

colnames(clinical_data)[grep("follow", colnames(clinical_data), ignore.case = TRUE)]

# Create survival dataframe
survival_df <- data.frame(
  patient_id = clinical_data$bcr_patient_barcode,
  vital_status = clinical_data$vital_status,
  days_to_death = clinical_data$days_to_death,
  days_to_last_followup = clinical_data$days_to_last_follow_up
)

# Check survival status distro
table(survival_df$vital_status)

# create time-to-event variable
survival_df$time <- ifelse(survival_df$vital_status == 'Dead',
                           survival_df$days_to_death,
                           survival_df$days_to_last_followup)

# Create event indicator
survival_df$event <- ifelse(survival_df$vital_status == "Dead",1,0)

# Remove patitnts with missing time data
survival_df <- survival_df[!is.na(survival_df$time), ]
print(paste("Patients for survival analysis: ", nrow(survival_df)))

# Create survival object
surv_object <- Surv(time = survival_df$time, event = survival_df$event)

#Fit Kaplan-Meier survival curve
km_fit <- survfit(surv_object ~ 1, data=survival_df)

# Plot overall curve with risk table
ggsurvplot(km_fit,
           data = survival_df,
           risk.table = TRUE,
           title = "Overall Survival - TCGA BRCA Patients",
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           legend = "none")

# Create enhanced survival plot
enhanced_plot <- ggsurvplot(km_fit,
                            data = survival_df,
                            risk.table = TRUE,
                            conf.int = TRUE,           # Add confidence interval
                            pval = TRUE,              # Add p-value (for overall survival)
                            title = "Overall Survival Analysis - TCGA BRCA Cohort",
                            xlab = "Time (Days)",
                            ylab = "Survival Probability",
                            legend = "none",
                            ggtheme = theme_minimal(), # Cleaner theme
                            risk.table.height = 0.25,  # Adjust risk table size
                            surv.median.line = "hv",   # Add median survival line
                            tables.theme = theme_cleantable())

# Add customizations
enhanced_plot$plot <- enhanced_plot$plot +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))

# Display the enhanced plot
print(enhanced_plot)

# Calculating the key survival statistics
summary_stats <- summary(km_fit)
median_survival <- summary(km_fit)$table['median']
print('survival summary statustics:')
print(paste('Median survival time:', round(median_survival, 2), 'days'))
print(paste('Number of patients:', km_fit$n))
print(paste("Number of events(death):", sum(survival_df$event)))

# Survival probabilities at key time points
time_points <- c(365, 739, 1825)
surv_at_times <- summary(km_fit, time=time_points)
print("Survival probabilities at key time points:")
print(surv_at_times$surv)

# Calculating the actual observed survival rates
observed_survival <- summary(km_fit)

# Extract the actual time points and survival probabilities
actual_data <- data.frame(
  time = observed_survival$time,
  survival = observed_survival$surv,
  patients_Remaining = observed_survival$n.risk
)
# Show the last 10 time points (most recent observations)
print("Last 10 observed time points:")
tail(actual_data, 10)

# Find the maximum observed time
max_observed_time <- max(observed_survival$time)

# Create survival plot limited to observed data only
observed_plot <- ggsurvplot(km_fit,
                            data = survival_df,
                            risk.table = TRUE,
                            conf.int = TRUE,
                            xlim = c(0, max_observed_time),  # Only show observed range
                            title = "Observed Survival (No Extrapolation) - TCGA BRCA",
                            xlab = "Time (Days)",
                            ylab = "Survival Probability",
                            legend = "none")

print(observed_plot)

print(paste("Maximum observed follow-up time:", max_observed_time, "days"))
print(paste("Survival probability at last observation:", 
            round(tail(observed_survival$surv, 1), 4)))

# Print the subtype information available
print("Available subtype columns:")
colnames(clinical_data)[grep("subtype|pam50|brca_subtype", colnames(clinical_data),
                             ignore.case = TRUE)]

# Check the pam50 subtype distro
if ("paper_BRCA_Subtype_PAM50" %in% colnames(clinical_data)) {
  print("PAM50 Subtype Distro:")
  table(clinical_data$paper_BRCA_Subtype_PAM50, useNA='always')
}

# Add subtype information to survival dataframe
survival_df$subtype <- clinical_data$paper_BRCA_Subtype_PAM50[match(
  survival_df$patient_id, clinical_data$bcr_patient_barcode)]
# Remove patients with missing subtype data
subtype_survival <- survival_df[!is.na(survival_df$subtype), ]
print("Patients with subtype information:")
table(subtype_survival$subtype)

# Create survival object for subtypes
subtype_surv_object <- Surv(time = subtype_survival$time,
                            event = subtype_survival$event)
# Fit Kaplan-Meier curves for each subtype
km_subtype_fit <- survfit(subtype_surv_object ~ subtype, data = subtype_survival)
# Plot subtype survival comparison
subtype_plot <- ggsurvplot(km_subtype_fit,
                           data = subtype_survival,
                           risk.table = TRUE,
                           conf.int = FALSE,
                           pval = TRUE,
                           title = "Survival by Breast Cancer Subtypes - TCGA BRCA",
                           xlab = "Time (Days)",
                           ylab = "Survival Probability",
                           legend.title = "PAM50 Subtype",
                           legend.labs = c("Basal", "HER2", "LumA", "LumB", "Normal"),
                           risk.table.height = 0.25)

print(subtype_plot)

# Get median survival for each subtypes
subtype_medians <- surv_median(km_subtype_fit)
print("Median survival by subtype:")
print(subtype_medians)
# Count events by subtype
print("Number of deaths by subtype:")
table(subtype_survival$subtype, subtype_survival$event)

# Check if certain subtypes have shorter follow-up times
subtype_survival %>%
  group_by(subtype) %>%
  summarize(
    max_followup = max(time),
    median_followup = median(time),
    n_patients = n()
  )

ggplot(subtype_survival, aes(x = subtype, y = time, color = as.factor(event))) +
  geom_boxplot() +
  labs(title = "Follow-up Time and Events by Subtype",
       x = "Subtype", 
       y = "Time (Days)",
       color = "Event (0=Alive, 1=Dead)") +
  theme_minimal()

# Find the maximum time where all subtypes have at least 10% of patients remaining
subtype_summary <- subtype_survival %>%
  group_by(subtype) %>%
  summarize(max_time = max(time),
            q10_time = quantile(time, 0.1))

common_followup <- min(subtype_summary$max_time)
print(paste("Common follow-up period for all subtypes:", common_followup, "days"))

# Create truncated dataset
truncated_survival <- subtype_survival[subtype_survival$time <= common_followup, ]

print("Patients remaining after truncation:")
table(truncated_survival$subtype)

# Create survival curves with balanced follow-up
trunc_surv_object <- Surv(time = truncated_survival$time, 
                          event = truncated_survival$event)

km_trunc_fit <- survfit(trunc_surv_object ~ subtype, data = truncated_survival)

# Plot truncated survival comparison
trunc_plot <- ggsurvplot(km_trunc_fit,
                         data = truncated_survival,
                         risk.table = TRUE,
                         pval = TRUE,
                         title = "Survival by Subtype (Balanced Follow-up) - TCGA BRCA",
                         xlab = "Time (Days)",
                         ylab = "Survival Probability",
                         legend.title = "PAM50 Subtype",
                         legend.labs = c("Basal", "HER2", "LumA", "LumB", "Normal"),
                         risk.table.height = 0.25)

print(trunc_plot)

# Get new median survivals
trunc_medians <- surv_median(km_trunc_fit)
print("Adjusted median survival by subtype:")
print(trunc_medians)

# Identify subtypes with sufficient long-term data (max follow-up > 4000 days)
adequate_subtypes <- subtype_summary[subtype_summary$max_time > 4000, ]$subtype
print("Subtypes with adequate long-term follow-up:")
print(adequate_subtypes)

# Analyze only these subtypes
adequate_survival <- subtype_survival[subtype_survival$subtype %in% adequate_subtypes, ]

# Create survival curves for adequate subtypes only
adequate_plot <- ggsurvplot(survfit(Surv(time, event) ~ subtype, data = adequate_survival),
                            data = adequate_survival,
                            risk.table = TRUE,
                            pval = TRUE,
                            title = "Survival - Subtypes with Adequate Follow-up",
                            xlab = "Time (Days)",
                            ylab = "Survival Probability")

print(adequate_plot)

# Verify the adjustment
print("Maximum follow-up in truncated data:")
max_followup_check <- truncated_survival %>%
  group_by(subtype) %>%
  summarize(max_time = max(time))

print(max_followup_check)

# Check if survival rankings changed
print("Original vs Adjusted Median Survival:")
comparison <- data.frame(
  Subtype = c("Basal", "Her2", "LumA", "LumB", "Normal"),
  Original_Median = c(7455, 6456, 3945, 3941, 4267),
  Adjusted_Median = c(3472, 3063, 3462, 3941, 4267)
)
print(comparison)

# Check if survival differences are statistically significant

# Log-rank test for subtype differences
subtype_test <- survdiff(Surv(time, event) ~ subtype, data = truncated_survival)
print("Statistical significance of subtype differences:")
print(subtype_test)

# Pairwise comparisons between key subtypes
basal_vs_her2 <- survdiff(Surv(time, event) ~ subtype, 
                          data = truncated_survival[truncated_survival$subtype %in% c("Basal", "Her2"), ])
print("Basal vs HER2 comparison:")
print(basal_vs_her2)

# Creating a comprehensive summary
final_summary <- data.frame(
  Subtype = c("Basal", "HER2", "LumA", "LumB", "Normal"),
  Patients = c(190, 81, 553, 208, 40),
  Deaths = c(28, 16, 64, 33, 7),
  Death_Rate = c(28/190, 16/81, 64/553, 33/208, 7/40),
  Median_Survival_Days = c(3472, 3063, 3462, 3941, 4267),
  Median_Survival_Years = round(c(3472, 3063, 3462, 3941, 4267)/365, 1)
)

print("Final Project Summary:")
print(final_summary)

print("Key Insights:")
print("1. Statistical tests show no strong evidence of survival differences between subtypes")
print("2. HER2 shows trend toward worse survival but needs larger sample")
print("3. Normal tissue samples show best survival (as expected)")
print("4. Proper follow-up adjustment was crucial to avoid wrong conclusions")
