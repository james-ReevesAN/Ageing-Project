#Load libraries
library(readr)
library(dplyr)
library(glmnet)
library(caret)
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
Agestr <- c("BATF"   ,  "CD160" ,   "CCL13"   , "FOXP3"   , "CD40LG" ,  "CD27"  ,   "BTLA"  ,   "IRAK1"  ,  "KLRF1" ,   "NLRP3"  ,  "TIGIT"  ,  "LEF1"  ,   "GZMB"   ,  "HLA-DPA1", "CX3CR1"  , "IRF1"   ,  "TNFRSF1B"  ,   "CXCR4" , "NT5E" , "CD163" , "NLRP3")


Agestr <- c("FOXP3" ,   "BATF"   ,  "CCL13"  ,  "CD160"  ,  "CD40LG"  , "KLRF1"   , "IRAK1" ,   "LEF1"  ,   "BTLA" ,    "HLA-DPA1" ,"IRF1")
#Agestr <- c("BACH2","TIGIT","KLRF1","PRF1","SESN3","LEF1","CCR7","B3GAT1","COX16","MYC", "GZMB" , "CD79A" , "FOXP3" , "CXCR4" , 'HLA-DPA1', "CX3CR1" , "CD27" , "CX3CR1" , "NT5E" , "CD163" , "TNFRSF1B" , "IL15", "NLRP3", "BTLA", "IFNG" , "IL7R" , "CD160",
#     "RORC", "CD97", "CD40LG" , "IRAK1", "ITGA6" , "IRF1" , "CCL13" , "BATF")

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

set.seed(123)

n <- length(survival_status)
K <- 5  # number of folds
folds <- createFolds(survival_status, k = K, list = TRUE)

predictions <- rep(NA, n)
true_labels <- survival_status

selected_gene_list <- vector("list", K)
selected_coef_list <- vector("list", K)

for (k in seq_len(K)) {
  # Train/test split
  test_idx <- folds[[k]]
  train_idx <- setdiff(seq_len(n), test_idx)
  
  x_train <- expr_mat[train_idx, , drop = FALSE]
  y_train <- survival_status[train_idx]
  x_test  <- expr_mat[test_idx, , drop = FALSE]
  
  # Ensure matrices
  x_train <- as.matrix(x_train)
  x_test  <- as.matrix(x_test)
  
  # Cross-validated LASSO within training set
  cv_model <- cv.glmnet(
    x_train, y_train,
    alpha = 1,
    family = "binomial",
    type.measure = "class",
    nfolds = max(2, length(train_idx) - 1),
    grouped = FALSE
  )
  
  best_lambda <- cv_model$lambda.min
  
  # Fit final model
  final_model <- glmnet(
    x_train, y_train,
    alpha = 1,
    lambda = best_lambda,
    family = "binomial"
  )
  
  # Prediction on test fold
  prob <- predict(final_model, s = best_lambda, newx = x_test, type = "response")
  predictions[test_idx] <- prob
  
  # Non-zero coefficients
  coef_k <- coef(final_model)
  nonzero_idx <- which(coef_k != 0)
  
  selected_gene_list[[k]] <- rownames(coef_k)[nonzero_idx]
  selected_coef_list[[k]] <- setNames(as.vector(coef_k[nonzero_idx]),
                                      rownames(coef_k)[nonzero_idx])
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


# Combine all coefficients into one data frame
coef_df <- do.call(rbind, lapply(seq_along(selected_coef_list), function(k) {
  coef_vec <- selected_coef_list[[k]]
  data.frame(
    Gene = names(coef_vec),
    Coefficient = as.numeric(coef_vec),
    Fold = k,
    stringsAsFactors = FALSE
  )
}))

# Drop intercept
coef_df <- coef_df %>% filter(Gene != "(Intercept)")

# Average coefficient per gene (without frequency)
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
  mutate(RiskGroup = ifelse(RiskScore > -5.1, "Low", "High" ))

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



