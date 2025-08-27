#Read in data
gene_data <- read_csv("~/R_Projects/NanoData/ANDonors/Normalized counts.nsolver.csv")
gene_data <- as.data.frame(gene_data)
clinical_data <- read_csv("~/R_Projects/NanoData/ANDonors/Groups-Donors copy.csv")
clinical_data <- as.data.frame(clinical_data)
ctrlGenes <- c("ERCC_00002.1","ERCC_00019.1","ERCC_00034.1","ERCC_00035.1","ERCC_00041.1","ERCC_00076.1","ERCC_00092.1","ERCC_00096.1","ERCC_00098.1","ERCC_00112.1","ERCC_00117.1",
               "ERCC_00126.1","ERCC_00144.1","ERCC_00154.1")
gene_data <- gene_data[!(rownames(gene_data) %in% ctrlGenes), ]
remove(ctrlGenes)

####Call Ageing genes in nSolver genenames
#Adapt to genes of interest:
Agestr <- c("BATF"   ,  "CD160" ,   "CCL13"   , "FOXP3"   , "CD40LG" ,  "CD27"  ,   "BTLA"  ,   "IRAK1"  ,  "KLRF1" ,   "NLRP3"  ,  "TIGIT"  ,  "LEF1"  ,   "GZMB"   ,  "HLA-DPA1", "CX3CR1"  , "IRF1"   ,  "RORC"  ,   "CXCR4" )
Age_df <- gene_data[gene_data$Name %in% Agestr, ]
remove(Agestr)

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

# Final expression matrix and survival vector
expr_mat <- expr_mat[ , !is.na(matched_indices)]  # Remove columns with unmatched IDs, if any
survival_status <- as.numeric(clinical_data$Survival)

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

#Step 4: Leave-One-Out Cross-Validation (LOOCV)
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
  
  # Cross-validated LASSO to find best lambda
  cv_model <- glmnet::cv.glmnet(x_train, y_train, alpha = 0, family = "binomial", type.measure = "class", nfolds = (nrow(x_train)-1), grouped = FALSE)
  best_lambda <- cv_model$lambda.min
  
  # Fit final model, save genes and coefficients and predict on test sample
  final_model <- glmnet(x_train, y_train, alpha = 0, lambda = best_lambda, family = "binomial")
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

# Plot ROC
plot(roc_obj, col = "blue", main = "LOOCV ROC Curve")

# Extract non-zero coefficients (excluding intercept)
coef_matrix <- coef(final_model)
non_zero_coefs <- coef_matrix[which(coef_matrix != 0), ]
non_zero_genes <- rownames(coef_matrix)[which(coef_matrix != 0)]

# Drop intercept
non_zero_genes <- non_zero_genes[-1]
non_zero_coefs <- non_zero_coefs[-1]

# Create a data frame of gene importance
gene_importance <- data.frame(
  Gene = non_zero_genes,
  Coefficient = as.numeric(non_zero_coefs)
)



# Rank by absolute coefficient
gene_importance <- gene_importance %>%
  arrange(desc(abs(Coefficient)))

#Step 6: Compute Patient Risk Scores
gene_Str <- gene_importance[,1]

# Get coefficients from final model (excluding intercept)
coef_vector <- coef(final_model)[-1, ]  # Remove intercept
gene_names_model <- rownames(coef(final_model))[-1]

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
library(ggplot2)

ggplot(patient_scores, aes(x = factor(Survival), y = RiskScore, fill = factor(Survival))) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(
    title = "Risk Scores by Survival Status",
    x = "Survival Status (0 = Survived, 1 =Died )",
    y = "Risk Score"
  ) +
  scale_fill_manual(values = c("green","red")) +
  theme_minimal()

risk_str <- as.vector(patient_scores[,2])

# AUC and ROC Curve for T.test model
#roc_obj <- multiclass.roc(risk_str, survival_status)
roc_obj <- roc(patient_scores$Survival, patient_scores$RiskScore)
auc_value <- auc(roc_obj)
cat("LOOCV AUC:", round(auc_value, 3), "\n")

# Plot ROC
plot(roc_obj, col = "blue", main = "LOOCV ROC Curve")


#Plot ROC curve based on Risk score
predicted_classes <- ifelse(predictions > 0.5, 1, 0)


# Accuracy
accuracy <- mean(predicted_classes == true_labels)
cat("LOOCV Accuracy:", round(accuracy, 3), "\n")

live <- patient_scores[patient_scores$Survival == 0,]
dead <- patient_scores[patient_scores$Survival == 1,]


NormCounts.PCA <- Age_df[2:ncol(Age_df)]
NormCounts.PCA <- t(NormCounts.PCA)
rownames(NormCounts.PCA) <- colnames(Age_df[2:ncol(Age_df)])


#EXTRA: Classic PCA- make the PCA from the data and organise it as needed
library(ggfortify)
library(ggplot2)

pca_nSolver <- prcomp(NormCounts.PCA[1:110,], scale. = TRUE)
autoplot(pca_nSolver, data = clinical_data, colour= 'Survival')
ClassicPrincComp <- pca_nSolver[["rotation"]]
write.csv(ClassicPrincComp, "/Users/jamesreeves/R_Projects/NanoData/ANDonors/ClassicPrincComp.nsolver.csv")

################################################################################## UMAP  ############################################################################################################################## 
clinical_data <- read_csv("~/R_Projects/NanoData/ANDonors/Groups-Donors copy.csv")
clinical_data <- as.data.frame(clinical_data)

clinical_data <- clinical_data[,1:9]


#Normnalise Data:
NormCounts.PCA[, 1: ncol(NormCounts.PCA)] <- scale(NormCounts.PCA[, 1:ncol(NormCounts.PCA)])


#Set seed for random values
set.seed(42)

#Run UMAP Algorithmn
NormCounts.PCA[,1:ncol(NormCounts.PCA)]  |>
  uwot::umap(
    n_neighbors = 6, n_threads = (parallel::detectCores()-2),
    min_dist = 0.025,
    n_sgd_threads = (parallel::detectCores()-2), 
    verbose = TRUE
  ) -> umap_exprn

umap_exprn <- data.table::as.data.table(umap_exprn)
colnames(umap_exprn) <- c("UMAP_X", "UMAP_Y")


#Join data from clinical & GeneCounts
NormCounts.PCA$SAMPLE <- colnames(Age_df[2:ncol(Age_df)])
join <- dplyr::inner_join(clinical_data, NormCounts.PCA)
#cbind umap values to exprn data.table
join <- cbind(join,umap_exprn)
remove(umap_exprn)

umap_del <- c("UMAP_X" , "UMAP_Y")
join <- join[,!(colnames(join) %in% umap_del) ]
###########UMAP Plots
# Convert 1/0 to Yes/No
join <- join %>%
  mutate(SurvivalFive = ifelse(SurvivalFive == 1, "Yes", "No"))
join$DonorAge <- as.numeric(join$DonorAge)

ggplot(join, aes(x = UMAP_X, y = UMAP_Y, color = SurvivalFive)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c( "Yes" = "green3", "No" = "orangered")) +
  theme_minimal()

ggplot(join, aes(x = UMAP_X, y = UMAP_Y, color = DonorAge)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient2(
    low = "#FFA46B",       # Low scores
    mid = "#F1EAE4",      # Midpoint (score = 0)
    high = "#5830F7",       # High scores
    midpoint = 35        # This sets 0 as the central color
  ) +
  theme_minimal() 

ggplot(join, aes(x = UMAP_X, y = UMAP_Y, color = Disease)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c( "MDS" = "#FFA46B", "AML" = "#F1EAE4", "ALL" = "#5830F7")) +
  theme_minimal()

################################################################################## SPECTRE::T-SNE  ############################################################################################################################## 
select_cols <- colnames(NormCounts.PCA)
rownames(NormCounts.PCA) <- colnames(Age_df[2:ncol(Age_df)])

#BROKEN
#Spectre::run.fitsne( NormCounts.PCA, use.cols = select_cols,
#  dims= 2,
#   theta = 0.5,
#    intervals_per_integer = 1,
#  ) -> tsne_exprn

tsne <- tsne::tsne(NormCounts.PCA)

tsne <- Rtsne::Rtsne(NormCounts.PCA)
tsne <- as.data.frame(tsne$Y)  
join <- cbind(join, tsne)

NormCounts.PCA <- as.data.frame(NormCounts.PCA)
colnames(NormCounts.PCA[,18:19]) <- c("X", "Y")

#Join data from 
NormCounts.PCA$SAMPLE <- rownames(NormCounts.PCA)
join <- dplyr::inner_join(clinical_data, NormCounts.PCA)
join$DonorAge <- as.numeric(join$DonorAge)

ggplot(join, aes(x = `1`, y = `2`, color = Disease)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c( "MDS" = "#FFA46B", "AML" = "#F1EAE4", "ALL" = "#5830F7")) +
  theme_minimal()

ggplot(join, aes(x = `1`, y = `2`, color = DonorAge)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient2(
    low = "#FFA46B",       # Low scores
    mid = "#F1EAE4",      # Midpoint (score = 0)
    high = "#5830F7",       # High scores
    midpoint = 35       # This sets 0 as the central color
  ) +
  theme_minimal() 

ggplot(join, aes(x = `1`, y = `2`, color = SurvivalFive)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c( "Yes" = "green3", "No" = "orangered")) +
  theme_minimal()
