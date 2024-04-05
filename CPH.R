# Load libraries
library('survival')
library(survex)
library(caret)
library(Hmisc)

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

################################################################################
# Define evaluation metrics

# Define function for calculating prediction error for 5-year RMST
calculate_pred_error <- function(coxphFit, data_test, time_cutoff = 1826.25) {
  # Calculate survival probabilities
  surv_prob <- survfit(coxphFit, newdata = data_test)
  surv <- surv_prob$surv
  
  # Calculate RMST
  rmst_i <- data.frame(matrix(ncol = 0, nrow = 0))
  
  # Define the time points for which to calculate the survival probabilities
  time_list <- c()
  
  for (i in surv_prob$time) {
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
    row_num <- which(surv_prob$time == t)
    s_t <- surv_prob$surv[row_num, ]
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

calc_ibs <- function(coxphFit, data_test) {
  # Extract the explanatory variables from the test dataset
  test_data_explanatory <- data_test[, !names(data_test) %in% c("time", "status")]
  
  # Predict survival probabilities for the test dataset
  cph_exp <- explain(coxphFit)
  surv <- cph_exp$predict_survival_function(coxphFit, test_data_explanatory, times = cph_exp$times)
  
  # Extract the actual survival times and censoring indicators from the test dataset
  y <- Surv(data_test$time, data_test$status)
  
  # Calculate the integrated Brier score to evaluate the model on the test set
  ibs <- integrated_brier_score(y, surv = surv, times = cph_exp$times)
  
  # Return the integrated Brier score
  return(ibs)
}

################################################################################
# Evaluation and Prediction

# Initialize variables to store results
cindex_results <- numeric(20)
ibs_results <- numeric(20)
prediction_error_results <- numeric(20)

# Repeat the process 20 times
for (i in 1:20) {
  # Make predictions on the test set using the Cox model
  pred_cox_test <- predict(coxphFit, newdata = data_test, type = "survival")
  
  #pred_cox_test <- survfit(coxphFit, newdata = data_test)$surv
  
  # Calculate c-index
  # Generate predicted survival probabilities
  predicted_survival <- predict(coxphFit, newdata = data_test, type = "survival")
  # Calculate the concordance index
  cindex_results[i] <- rcorr.cens(predicted_survival, Surv(data_test$time, data_test$status))[1]
  
  # Calculate IBS
  ibs <- calc_ibs(coxphFit, data_test)
  ibs_results[i] <- ibs
  
  # Calculate prediction error
  pred_error <- calculate_pred_error(coxphFit, data_test)[[2]]
  prediction_error_results[i] <- pred_error
}

# Calculate mean c-index and mean IBS
mean_cindex <- mean(cindex_results)
mean_ibs <- mean(ibs_results)
mean_pred_error <- mean(prediction_error_results)

# Report the results
cat("Mean c-index:", round(mean_cindex, 3), "\n")
cat("Mean IBS:", round(mean_ibs, 3), "\n")
cat("Mean prediction error:", round(mean_pred_error, 3), "\n")

# Plot boxplot of c-index
boxplot(cindex_results, main = "Distribution of C-index over 20 Replicates")

# Plot boxplot of Integrated Brier score
boxplot(ibs_results, main = "Distribution of Integrated Brier Score over 20 Replicates")

# Plot boxplot of prediction error
boxplot(prediction_error_results, main = "Distribution of Prediction Error over 20 Replicates")


# Combine the lists into one string
combined_data <- paste("# CPH\n",
                       "cph_c_index=[", paste(cindex_results, collapse = ", "), "]\n",
                       "cph_ibs=[", paste(ibs_results, collapse = ", "), "]\n",
                       "cph_pe=[", paste(prediction_error_results, collapse = ", "),"]", sep = "")

# Save the combined data to a text file
writeLines(combined_data, "CPH_scores.txt")

