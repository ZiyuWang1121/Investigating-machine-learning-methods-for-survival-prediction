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
#brca <- read.csv("brca_filtered.csv")

# Split the data into train/test sets
#train <- sample(1:nrow(brca), round(nrow(brca) * 0.80))

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

################################################################################
# Parameter tunning

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
    cindex <- get.cindex(model$yvar[,1], model$yvar[,2], model$predicted.oob)
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

################################################################################
# Define evaluation metrics

# Define function for calculating prediction error for 5-year RMST
calculate_pred_error <- function(final_model, data_test, time_cutoff = 1826.25) {
  # Calculate survival probabilities
  pred <- predict(final_model, data_test)
  surv <- pred$survival
  
  # Calculate RMST
  rmst_i <- data.frame(matrix(ncol = 0, nrow = 0))
  
  time_list <- c()
  
  for (i in pred$time) {
    if (i <= (time_cutoff + 1)) {
      time_list <- c(time_list, i)
    }
  }
  
  if (!(time_cutoff %in% time_list)) {
    time_list <- c(time_list, time_cutoff)
  }
  
  for (i in 1:(length(time_list) - 1)) {
    t <- time_list[i]
    t_plus_1 <- time_list[i + 1]
    row_num <- which(pred$time == t)
    s_t <- pred$survival[, row_num]
    mu <- (t_plus_1 - t) * s_t
    
    # Append mu to rmst_i
    rmst_i <- rbind(rmst_i, mu)
  }
  
  rmst <- data.frame(matrix(0, nrow = 1, ncol = ncol(rmst_i)))
  colnames(rmst) <- colnames(rmst_i)
  for (i in 1:ncol(rmst_i)) {
    rmst[1, i] <- sum(rmst_i[, i])
  }
  
  rmst.T <- t(rmst)
  
  test <- data_test[, c(1, 2)]
  row.names(test) <- NULL
  test$rmst <- rmst.T
  test <- test[test$time <= time_cutoff, ]
  
  train <- data_train[, c(1, 2)]
  train <- train[train$time <= time_cutoff, ]
  
  # Kaplan-Meier estimator of the censoring distribution
  # fitted to training dataset
  kmf <- survfit(Surv(train$time, 1 - train$status) ~ 1)
  
  cen <- data.frame(time = kmf$time, surv_prob = kmf$surv)
  cen <- cen[cen$time <= time_cutoff, ]
  
  # Calculate the time intervals col
  cen$interval <- cbind(cen$time, c(cen$time[-1], time_cutoff))
  
  survival_probabilities <- c()
  
  # Iterate over each time point
  for (i in test$time) {
    # Find the corresponding interval in the 'interval' column of 'cen'
    for (j in 1:nrow(cen)) {
      if (i >= cen$interval[j, 1] && i < cen$interval[j, 2]) {
        # Append the survival probability to the list
        survival_probabilities <- c(survival_probabilities, cen$surv_prob[j])
        break
      }
    }
  }
  
  # Add the list of survival probabilities as a new column to the DataFrame
  test$KM_estimate <- survival_probabilities
  
  up_sum <- 0
  low_sum <- 0
  
  for (i in 1:nrow(test)) {
    # 5-year RMST prediction error
    up <- (1 / test$KM_estimate[i]) * test$status[i] * abs(test$time[i] - test$rmst[i])
    low <- (1 / test$KM_estimate[i]) * test$status[i]
    up_sum <- up_sum + up
    low_sum <- low_sum + low
  }
  
  pred_error <- up_sum / low_sum
  
  return(list(rmst, pred_error))
}

calc_ibs <- function(final_model, data_test) {
  # Extract the actual survival times and censoring indicators from the test dataset
  y <- Surv(data_test$time, data_test$status)
  
  # Calculate the integrated Brier score to evaluate the model on the test set
  ibs <- integrated_brier_score(y, surv = cbind(pred_test$survival, pred_test$survival.oob), times = pred_test$time.interest)
  
  # Return the integrated Brier score
  return(ibs)
}

################################################################################
# Evaluation and Prediction

# Initialize vectors to store results
cindex_results <- numeric(20)
#brier_scores <- numeric(20)
integrated_brier_scores <- numeric(20)
prediction_error_results <- numeric(20)

# mtry = 103, nodesize = 20
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
  cindex_test <- get.cindex(final_model$yvar[,1], final_model$yvar[,2], final_model$predicted.oob)
  cindex_results[i] <- cindex_test
  
  # Calculate Brier score by Kaplan-Meier method
  #bs_km <- get.brier.survival(final_model, cens.model = "km")$brier.score
  #brier_scores[i] <- bs_km
  
  # Calculate IBS
  ibs <- calc_ibs(final_model, data_test)
  integrated_brier_scores[i] <- ibs
  
  # Calculate prediction error
  pred_error <- calculate_pred_error(final_model, data_test)[[2]]
  prediction_error_results[i] <- pred_error
}

# Calculate mean c-index
mean_cindex <- mean(cindex_results)
mean_cindex

# Calculate mean integrated Brier score
mean_ibs <- mean(integrated_brier_scores)
mean_ibs

mean_pred_error <- mean(prediction_error_results)
mean_pred_error

# Plot boxplot of c-index
boxplot(cindex_results, main = "Distribution of C-index over 20 Replicates")

# Plot boxplot of Brier score
#boxplot(brier_scores, main = "Distribution of Brier Score over 20 Replicates")

# Plot boxplot of Integrated Brier score
boxplot(integrated_brier_scores, main = "Distribution of Integrated Brier Score")

# Plot boxplot of prediction error
boxplot(prediction_error_results, main = "Distribution of Prediction Error over 20 Replicates")

# Combine the lists into one string
combined_data <- paste("# RSF\n",
                       "rsf_c_index=[", paste(cindex_results, collapse = ", "), "]\n",
                       "rsf_ibs=[", paste(integrated_brier_scores, collapse = ", "), "]\n",
                       "rsf_pe=[", paste(prediction_error_results, collapse = ", "),"]", sep = "")

# Save the combined data to a text file
writeLines(combined_data, "RSF_scores.txt")
