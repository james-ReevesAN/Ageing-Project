#Load libraries
library(readr)
library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(pheatmap)
library(survival)
library(survminer)


#Read in data
gene_data <- read_csv("~/R_Projects/NanoData/ANDonors/Normalized counts.nsolver.csv")
gene_data <- as.data.frame(gene_data)
clinical_data <- read_csv("~/R_Projects/NanoData/ANDonors/Groups-Donors copy.csv")
clinical_data <- as.data.frame(clinical_data)

#Adapt to genes of interest:
Agestr <- c("BATF"   ,  "CD160" ,   "CCL13"   , "FOXP3"   , "CD40LG" ,  "CD27"  ,   "BTLA"  ,   "IRAK1"  ,  "KLRF1" ,   "NLRP3"  ,  "TIGIT"  ,  "LEF1"  ,   "GZMB"   ,  "HLA-DPA1", "CX3CR1"  , "IRF1"   ,  "RORC"  ,   "CXCR4" )

Age_df <- gene_data[gene_data$Name %in% Agestr, ]
remove(Agestr)

#####################
############# ADD PMM USING MICE HERE 
#####################

# Step 2: Prepare gene expression matrix 
# Extract gene names (rows = genes, columns = patients)
gene_names <- Age_df$Name[1:nrow(Age_df)]
expression_matrix <- Age_df[1:nrow(Age_df), 4:ncol(Age_df)]  # Assuming column C onward are patient IDs

# Convert to matrix and set row and column names
expr_mat <- as.matrix(expression_matrix)
rownames(expr_mat) <- gene_names

# Transpose: rows = patients, columns = genes
expr_mat <- t(expr_mat)

#Step 3: Prepare clinical outcome data
# Assuming patient IDs are in column A and binary survival status in column G
clinical_data <- clinical_data [1:110,] %>%
  select(PatientID = 1, Survival = 8) 


# Match using last 8 characters of patient/sample IDs
# Find matching patient IDs
expr_ids <- rownames(expr_mat)
clinical_ids <- clinical_data$PatientID

# Reorder clinical data to match expression data column order
matched_indices <- match(expr_ids, clinical_ids)
clinical_data <- clinical_data[matched_indices, ]

# Remove unmatched rows (just in case)
expr_mat <- expr_mat[!is.na(matched_indices), ]
survival_status <- as.numeric(clinical_data$Survival[!is.na(matched_indices)])

# Step 3.5: Handle missing values
# Remove genes (columns) with any NA
na_cols <- apply(expr_mat, 2, function(x) any(is.na(x)))
expr_mat <- expr_mat[, !na_cols]

# Optional: remove samples with any NA (in case any still exist)
na_rows <- apply(expr_mat, 1, function(x) any(is.na(x)))
expr_mat <- expr_mat[!na_rows, ]
survival_status <- survival_status[!na_rows]
remove(na_cols)
remove(na_rows)
remove(clinical_ids)
remove(expr_ids)


#Step 4: K FOld Cross Validation Lasso 
predictions <- rep(NA, length(survival_status))
true_labels <- survival_status

#Create the lists for getting the list of the genes in the model
selected_gene_list <- list()
selected_coef_list <- list()

set.seed(123)  # For reproducibility

for (i in 1:length(survival_status)) {
  # Split into training and test
  x_train <- expr_mat[-i, , drop = FALSE]
  y_train <- survival_status[-i]
  x_test <- expr_mat[i, , drop = FALSE]
  
  # Ensure x_train is a matrix and has ≥2 columns
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  
  # K Fold Cross-validated LASSO to find best lambda
  cv_model <- glmnet::cv.glmnet(x_train, y_train, alpha = 1, family = "binomial", type.measure = "class", nfolds = (nrow(x_train)-1), grouped = FALSE)
  best_lambda <- cv_model$lambda.min
  
  # Fit final model, save genes and coefficients and predict on test sample
  final_model <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda, family = "binomial")
  prob <- predict(final_model, s = best_lambda, newx = x_test, type = "response")
  predictions[i] <- prob
  
  # Get non-zero coefficients
  coef_i <- coef(final_model)
  selected_genes <- rownames(coef_i)[which(coef_i != 0)]
  selected_coefs <- coef_i[which(coef_i != 0)]
  
  # Save to a list
  selected_gene_list[[i]] <- selected_genes
  selected_coef_list[[i]] <- selected_coefs
}


#Step 5: Evaluate Performance

# Binary predictions (threshold 0.5)
predicted_classes <- ifelse(predictions > 0.5, 1, 0)

# Accuracy
accuracy <- mean(predicted_classes == true_labels)
cat("LOOCV Accuracy:", round(accuracy, 3), "\n")

# AUC and ROC Curve
roc_obj <- roc(true_labels, predictions)
auc_value <- auc(roc_obj)
cat("LOOCV AUC:", round(auc_value, 3), "\n")


# Combine coefficients from all folds into one data frame
coef_df <- do.call(rbind, lapply(seq_along(selected_coef_list), function(i) {
  coef_vec <- selected_coef_list[[i]]
  genes <- selected_gene_list[[i]]
  
  # Make sure lengths match (and drop intercept if present)
  keep_idx <- genes != "(Intercept)"
  genes <- genes[keep_idx]
  coefs <- as.numeric(coef_vec[keep_idx])
  
  if (length(genes) == 0) return(NULL)
  
  data.frame(
    Gene = genes,
    Coefficient = coefs,
    Fold = i,
    stringsAsFactors = FALSE
  )
}))

# Summarise across folds
gene_importance <- coef_df %>%
  group_by(Gene) %>%
  summarise(
    Coefficient = mean(Coefficient, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(Coefficient)))
gene_importance <- as.data.frame(gene_importance)

# Rank by absolute coefficient
gene_importance <- gene_importance %>%
  arrange(desc(abs(Coefficient)))

#Step 6: Compute Patient Risk Scores
gene_Str <- gene_importance[,1]

# Get coefficients from final model (excluding intercept)
coef_vector <-  setNames(
  gene_importance$Coefficient,
  gene_importance$Gene
)
gene_names_model <- gene_importance$Gene

# Keep only genes with non-zero coefficients
selected_genes <- gene_names_model[coef_vector != 0]
selected_coefs <- coef_vector[coef_vector != 0]

# Subset expression matrix to only selected genes
expr_selected <- expr_mat[, selected_genes]

# Ensure dimensions align (some rows might need coercion)
expr_selected <- as.matrix(expr_selected)

# Compute score: matrix multiplication of gene expression and coefficients
risk_scores <- expr_selected %*% selected_coefs

# Attach patient IDs for reference
patient_scores <- data.frame(
  PatientID = rownames(expr_mat),
  RiskScore = as.numeric(risk_scores)
)

# Attach survival status to patient scores
patient_scores$Survival <- true_labels  # already ordered to match expr_mat

# Plot: Risk Score by Survival Status


ggplot(patient_scores, aes(x = factor(Survival), y = RiskScore, fill = factor(Survival))) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(
    title = "Risk Scores by Survival Status",
    x = "Survival Status (0 = Survived, 1 =Died )",
    y = "Risk Score"
  ) +
  scale_fill_manual(values = c("#00D600","#184918")) +
  theme_minimal()

# AUC and ROC Curve for T.test model
#roc_obj <- multiclass.roc(risk_str, survival_status)
roc_obj <- roc(patient_scores$Survival, patient_scores$RiskScore)
auc_value <- auc(roc_obj)
cat("LOOCV AUC:", round(auc_value, 3), "\n")

# Plot ROC
plot(roc_obj, col = "blue", main = "LOOCV ROC Curve")


# UMAP:
clinical_data <- read_csv("~/R_Projects/NanoData/ANDonors/Groups-Donors copy.csv")
clinical_data <- as.data.frame(clinical_data)
clinical_data <- clinical_data[,1:10]


NormCounts <- as.data.frame(expr_mat)
# Convert to numeric and log2 transform (+ pseudocount to avoid log2(0))
NormCounts.PCA <- as.data.frame(sapply(NormCounts, function(x) log2(as.numeric(x) + 1)))
# Restore row names (patient IDs)
rownames(NormCounts.PCA) <- rownames(NormCounts)
remove(NormCounts)
#Join data from clinical & GeneCounts
NormCounts.PCA$SAMPLE <- colnames(Age_df[4:ncol(Age_df)])
join <- dplyr::inner_join(clinical_data, NormCounts.PCA)

# Add Patient Risk scores
patient_scores <- patient_scores[,1:2]
colnames(patient_scores) <- c("SAMPLE" , "RiskScore")
join <- left_join(join, patient_scores)
#Assign Risk Group - High/Low
join <- join %>% 
  mutate(RiskGroup = ifelse(RiskScore > -5, "Low", "High" ))

# Ensure correct data types
join$SurvivalFive <- as.numeric(ifelse(join$SurvivalFive == 1, 0, 1))       # 1 = event, 0 = censored
join$SurvivalTime <- as.numeric(join$SurvivalTime) * 12
join$RiskGroup <- as.factor(join$RiskGroup)          # column with "High" / "Low"

# Build survival object
surv_obj <- Surv(time = join$SurvivalTime, event = join$SurvivalFive)

# Fit Kaplan-Meier model grouped by High/Low
km_fit_group <- survfit(surv_obj ~ RiskGroup, data = join)

# Plot group-wise Kaplan-Meier curves
survminer::ggsurvplot(
  km_fit_group,
  data = join,
  fun = "pct",
  risk.table.title = "Survival Rate:",
  palette = c("#184918","#00D600"),
  risk.table = TRUE,       # number at risk table
  pval = TRUE,             # log-rank test
  conf.int = TRUE,         # confidence intervals
  legend.title = "Risk Group:",
  legend.labs = c("High", "Low"),
  xlab = "Time in Months",
  ylab = "Survival Probability (%)"
)



