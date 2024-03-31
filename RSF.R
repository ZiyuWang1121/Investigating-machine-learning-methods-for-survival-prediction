# Load libraries
library('survival')
library('randomForestSRC')
library(ggplot2)
library(ggRandomForests)
library(caret)
library(recipes)
library(survex)

# set random state
set.seed(0)

# Load data
brca <- read.csv("brca.csv")

# Split the data into train/test sets
train <- sample(1:nrow(brca), round(nrow(brca) * 0.80))
data_train_orign <- brca[train, ]
data_test_orign <- brca[-train, ]

pca <- prcomp(data_train_orign[, -c(1,2)])

# Calculate the cumulative proportion of variance explained
cumulative_variance <- cumsum(pca$sdev^2) / sum(pca$sdev^2)

# Find the number of components required to explain 99% of the variance
num_components <- which(cumulative_variance >= 0.99)[1]

# Plot cumulative variance explained
ggplot(data.frame(Component = 1:length(cumulative_variance), 
                  Variance = cumulative_variance), 
       aes(x = Component, y = Variance)) +
  geom_line() +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "red") +
  labs(title = "Cumulative Variance Explained by Principal Components",
       x = "Number of Components",
       y = "Cumulative Variance Explained") +
  theme_minimal()

# Retain only the required number of components
#pca_train_data <- pca$x[, 1:num_components]

# Print the number of components required
print(paste("Number of components to explain 99% of variance:", num_components))

#biplot(pca)

# Summary of PCA
#summary(pca)

# Plot the variance explained by each principal component
#plot(pca)

# Apply PCA transformation to training and test sets
data_train_pca <- predict(pca , data_train_orign[, -c(1,2)])
data_test_pca <- predict(pca, data_test_orign[, -c(1,2)])

data_train = cbind(data_train_orign[, c(1, 2)], data_train_pca[, 1:num_components])
data_test = cbind(data_test_orign[, c(1, 2)], data_test_pca[, 1:num_components])


# Initialize variables to store results
cindex_scores <- c() 
best_cindex <- 0
best_params <- NULL

# Perform 5 repeats of 5-fold cross validation for parameter tuning
for (i in 1:5) {
  # Create 5-fold cross-validation folds
  folds <- createFolds(data_train$status, k = 5, returnTrain = TRUE)
  
  # Initialize variables to store fold results
  cindex_fold_scores <- c() 
  
  for (j in 1:5) {
    # Get training and validation data for the fold
    fold_train <- data_train[folds[[j]], ]
    fold_valid <- data_train[-folds[[j]], ]
    
    # Tune model on fold training data
    tuned_model <- tune.rfsrc(Surv(time, status) ~ .,
                              data = fold_train,
                              mtryStart = ncol(fold_train) / 2,
                              nodesizeTry = c(1:9, seq(10, 100, by = 5)),
                              ntreeTry = 100)
    # Refit model to fold training data using optimal parameters
    model <- rfsrc(Surv(time, status) ~ .,
                   data = fold_train,
                   mtry = tuned_model$optimal[2],
                   nodesize = tuned_model$optimal[1])
    
    # Evaluate model on fold validation data: Concordance Index
    pred <- predict(model, fold_valid)
    cindex <- get.cindex(fold_valid$time, fold_valid$status, -pred$predicted)
    cindex_fold_scores <- c(cindex_fold_scores, cindex)
    
    # Check if current C-index is the best so far
    if (cindex > best_cindex) {
      best_cindex <- cindex
      best_params <- list(nodesize = tuned_model$optimal[1], mtry = tuned_model$optimal[2])
    }
  }
  
  # Store average fold results
  cindex_scores <- c(cindex_scores, mean(cindex_fold_scores))
}

# Report the best c-index and corresponding parameters
cat("Best c-index:", best_cindex, "\n")
cat("Corresponding parameters: mtry =", best_params$mtry, ", nodesize =", best_params$nodesize, "\n")


# Initialize vectors to store results
cindex_results <- numeric(20)
#brier_scores <- numeric(20)
integrated_brier_scores <- numeric(20)

for (i in 1:20) {
  # Refit model to all training data using optimal parameters
  final_model <- rfsrc(Surv(time, status) ~ .,
                       data = data_train,
                       mtry = best_params$mtry,
                       nodesize = best_params$nodesize)
  
  # Make prediction on train set and c-index for evaluation
  #pred_train <- predict(final_model, data_train)
  #cindex_train <- get.cindex(data_train$time, data_train$status, -pred_train$predicted)
  
  # Make prediction on test set and c-index for evaluation
  pred_test <- predict(final_model, data_test)
  cindex_test <- get.cindex(data_test$time, data_test$status, -pred_test$predicted)
  cindex_results[i] <- cindex_test
  
  # Calculate Brier score by Kaplan-Meier method
  #bs_km <- get.brier.survival(final_model, cens.model = "km")$brier.score
  #brier_scores[i] <- bs_km
  
  # Calculate integrated Brier score
  model_exp <- explain(final_model)
  y <- model_exp$y
  times <- model_exp$times
  surv <- model_exp$predict_survival_function(final_model,  model_exp$data, times)
  ibs <- integrated_brier_score(y, surv = surv, times = times)
  integrated_brier_scores[i] <- ibs
}

# Calculate mean c-index
mean_cindex <- mean(cindex_results)
mean_cindex

# Calculate mean integrated Brier score
mean_ibs <- mean(integrated_brier_scores)
mean_ibs


# Plot boxplot of c-index
boxplot(cindex_results, main = "Distribution of C-index over 20 Replicates")

# Plot boxplot of Brier score
#boxplot(brier_scores, main = "Distribution of Brier Scoreover 20 Replicates")

# Plot boxplot of Integrated Brier score
boxplot(integrated_brier_scores, main = "Distribution of Integrated Brier Score")

###############################################################################

# Fit Cox proportional hazards model for baseline comparison
coxphFit <- coxph(Surv(time, status) ~ ., data = data_train,x=TRUE)

# Summary of Cox proportional hazards model
summary(coxphFit)

pred_cox_train <- predict(coxphFit, data_train)
# Calculate C-index for the Cox proportional hazards model on the training set
cindex_cox_train <- get.cindex(data_train$time, data_train$status, -pred_cox_train)
round(cindex_cox_train,3)

# Make predictions on the test set using the Cox model
pred_cox_test <- predict(coxphFit, newdata = data_test, type = "lp")
cindex_cox_test <- get.cindex(data_test$time, data_test$status, -pred_cox_test)
round(cindex_cox_test,3)

ibs_cox <- integrated_brier_score(data_test$status, surv = pred_cox_test, times = data_test$time)

