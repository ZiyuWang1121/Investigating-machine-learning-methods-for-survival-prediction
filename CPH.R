# Load libraries
library('survival')
library(survex)

# set random state
set.seed(0)

# Load data
brca <- read.csv("brca.csv")
#brca <- read.csv("brca_filtered.csv")

# 80% for training
# By default, createDataPartition does a stratified random split of the data.
train <- createDataPartition(brca$status, p = 0.8, list = FALSE)
data_train_orign <- brca[train, ]
data_test_orign <- brca[-train, ]

pca <- prcomp(data_train_orign[, -c(1,2)])

# Calculate the cumulative proportion of variance explained
cumulative_variance <- cumsum(pca$sdev^2) / sum(pca$sdev^2)

# Find the number of components required to explain 99% of the variance
num_components <- which(cumulative_variance >= 0.99)[1]

# Print the number of components required
print(paste("Number of components to explain 99% of variance:", num_components))

# Apply PCA transformation to training and test sets
data_train_pca <- predict(pca , data_train_orign[, -c(1,2)])
data_test_pca <- predict(pca, data_test_orign[, -c(1,2)])

data_train = cbind(data_train_orign[, c(1, 2)], data_train_pca[, 1:num_components])
data_test = cbind(data_test_orign[, c(1, 2)], data_test_pca[, 1:num_components])


# Fit Cox proportional hazards model for baseline comparison
coxphFit <- coxph(Surv(time, status) ~ ., data = data_train, model = TRUE, x = TRUE, y = TRUE)

# Summary of Cox proportional hazards model
summary(coxphFit)


# Initialize variables to store results
cindex_results <- numeric(20)
ibs_results <- numeric(20)

# Repeat the process 20 times
for (i in 1:20) {
  # Make predictions on the test set using the Cox model
  pred_cox_test <- predict(coxphFit, newdata = data_test, type = "lp")
  
  # Calculate c-index
  cindex_results[i] <- get.cindex(data_test$time, data_test$status, -pred_cox_test)
  
  model_exp <- explain(coxphFit)
  y <- model_exp$y
  times <- model_exp$times
  surv <- model_exp$predict_survival_function(coxphFit,  model_exp$data, times)
  ibs <- integrated_brier_score(y, surv = surv, times = times)
  ibs_results[i] <- ibs
  
}

# Calculate mean c-index and mean IBS
mean_cindex <- mean(cindex_results)
mean_ibs <- mean(ibs_results)

# Report the results
cat("Mean c-index:", round(mean_cindex, 3), "\n")
cat("Mean IBS:", round(mean_ibs, 3), "\n")

# Plot boxplot of c-index
boxplot(cindex_results, main = "Distribution of C-index over 20 Replicates")

# Plot boxplot of Integrated Brier score
boxplot(ibs_results, main = "Distribution of Integrated Brier Score")
