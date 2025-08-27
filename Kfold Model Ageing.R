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
Agestr <- c("BACH2","TIGIT","KLRF1","PRF1","SESN3","LEF1","CCR7","B3GAT1","COX16","MYC", "GZMB" , "CD79A" , "FOXP3" , "CXCR4" , 'HLA-DPA1', "CX3CR1" , "CD27" , "CX3CR1" , "NT5E" , "CD163" , "TNFRSF1B" , "IL15", "NLRP3", "BTLA", "IFNG" , "IL7R" , "CD160")

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




library(glmnet)
X <- expr_mat
y <- survival_status
set.seed(123)
run_lasso_cv <- function(X, y, k_folds = 10) {
  # Convert to proper formats
  X <- as.matrix(X)
  if (!is.numeric(y)) {
    y <- as.numeric(as.factor(y)) - 1  # Convert factor to 0/1
  }
  
  # Fit lasso with CV
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 1,
    nfolds = k_folds,
    type.measure = "class"
  )
  
  # Extract non-zero coefficients at best lambda
  coef_vec <- coef(cv_fit, s = "lambda.min")
  coef_vals <- as.vector(coef_vec)
  coef_names <- rownames(coef_vec)
  
  nonzero_mask <- coef_vals != 0
  nonzero_names <- coef_names[nonzero_mask]
  nonzero_values <- coef_vals[nonzero_mask]
  
  # Skip intercept (if present)
  if ("(Intercept)" %in% nonzero_names) {
    skip <- which(nonzero_names == "(Intercept)")
    nonzero_names <- nonzero_names[-skip]
    nonzero_values <- nonzero_values[-skip]
  }
  
  # Create dataframe
  result_df <- data.frame(
    Gene = nonzero_names,
    Coefficient = round(nonzero_values, 10),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}


genes <-  run_lasso_cv(X, y, 20)
gene_importance <- genes
genenames <- strsplit(genes, "=")
genenames <- strsplit(genes, ";")
