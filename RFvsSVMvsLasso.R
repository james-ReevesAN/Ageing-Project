# Load required libraries
library(randomForest)
library(tidyverse)
library(caret)
library(pROC)

expr_df <- as.data.frame(expr_mat)
expr_df$PatientID <- rownames(expr_df)
expr_df <- left_join(expr_df,clinical_data)

# Convert survival column to factor
expr_df$Survival <- as.factor(expr_df$Survival)
expr_df$Survival <- factor(expr_df$Survival, levels = c(0, 1), labels = c("No", "Yes"))
colnames(expr_df)[17] <- "HLADPA1"
# Split data into features and target
features <- expr_df[1:29]
target <- expr_df$Survival

# Train-test split
set.seed(123)
split_index <- createDataPartition(target, p = 0.8, list = FALSE)
train_x <- features[split_index, ]
train_y <- target[split_index]
test_x  <- features[-split_index, ]
test_y  <- target[-split_index]

# Normalize gene expression (Z-score scaling)
pre_proc <- preProcess(train_x, method = c("center", "scale"))
train_x <- predict(pre_proc, train_x)
test_x  <- predict(pre_proc, test_x)

# Combine training features and labels for caret training
dtrain <- data.frame(train_x, Survival = train_y)

tctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

# Train Random Forest model
set.seed(456)
rf_cv_model <- train(
  Survival ~ .,
  data = dtrain,
  method = "rf",
  trControl = tctrl,
  metric = "ROC",
  tuneLength = 5
)

# Train Logistic Regression model
set.seed(456)
glm_model <- train(
  Survival ~ .,
  data = dtrain,
  method = "glm",
  family = "binomial",
  trControl = tctrl,
  metric = "ROC"
)

# Train Support Vector Machine (radial) model
set.seed(456)
svm_model <- train(
  Survival ~ .,
  data = dtrain,
  method = "svmRadial",
  trControl = tctrl,
  metric = "ROC",
  tuneLength = 5
)

# Print performance
print(rf_cv_model)
print(glm_model)
print(svm_model)

# Predict probabilities on test data
rf_probs <- predict(rf_cv_model, newdata = test_x, type = "prob")
glm_probs <- predict(glm_model, newdata = test_x, type = "prob")
svm_probs <- predict(svm_model, newdata = test_x, type = "prob")

# Compute ROC curves
roc_rf  <- roc(test_y, rf_probs[, 2])
roc_glm <- roc(test_y, glm_probs[, 2])
roc_svm <- roc(test_y, svm_probs[, 2])

# Plot ROC curves
plot(roc_rf, col = "#5830F7", main = "Model Comparison")
lines(roc_glm, col = "#FFA46B")
lines(roc_svm, col = "#00D600")

legend("bottomright", legend = c("Random Forest", "Logistic Regression", "SVM"),
       col = c("#5830F7", "#FFA46B", "#00D600"), lwd = 0.75)

# Print AUCs
cat("Random Forest AUC:", auc(roc_rf), "\n")
cat("Logistic Regression AUC:", auc(roc_glm), "\n")
cat("SVM AUC:", auc(roc_svm), "\n")

# Variable importance plot
varImpPlot(rf_cv_model$finalModel, main = "Variable Importance - Top Genes")

# Create data frame for ggplot
roc_df <- bind_rows(
  data.frame(model = "Random Forest",  tpr = roc_rf$sensitivities, fpr = 1 - roc_rf$specificities),
  data.frame(model = "Logistic Regression", tpr = roc_glm$sensitivities, fpr = 1 - roc_glm$specificities),
  data.frame(model = "SVM", tpr = roc_svm$sensitivities, fpr = 1 - roc_svm$specificities)
)

# Plot with ggplot2
roc_df %>%
  ggplot(aes(x = fpr, y = tpr, color = model)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "ROC Curves - Model Comparison",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Model"
  )
