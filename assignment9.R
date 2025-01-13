# load useful R packages
library(ggplot2)
library(dplyr)
library(car)
library(stats)
# read in the data
load("assignment9_Selumetinib (AZD6244)_FLT3-ITD.RData")

data <- merge(clinical_data, drug_data, by = "sample")

AUC <- data$auc
FLT3 <- data$'FLT3-ITD'
#1.(i)
# Visualize the distributions
par(mfrow=c(1,2))
hist(AUC[which(FLT3==0)],main="FLT3 negative", xlab="AUC", ylab="Frequency", breaks = 10)
hist(AUC[which(FLT3==1)],main="FLT3 positive", xlab="AUC", ylab="Frequency", breaks = 10)

# Create boxplots
par(mfrow=c(1,1))
boxplot(AUC~FLT3, data=data, main="FLT3-ITD vs AUC", xlab="FLT3-ITD", ylab="AUC")

# Compare the statistics
library(skimr)
skim(AUC[which(FLT3==0)])
skim(AUC[which(FLT3==1)])

#1.(ii) Compare the two groups using a t-test
#Check assumptions
shapiro.test(AUC[which(FLT3==0)])
shapiro.test(AUC[which(FLT3==1)])
par(mfrow=c(1,2))
qqPlot(AUC[which(FLT3==0)], main="FLT3 negative")
qqPlot(AUC[which(FLT3==1)], main="FLT3 positive")
#Compare two groups using a t-test
t.test(AUC[which(FLT3==0)], AUC[which(FLT3==1)])
#Compare the two groups using a Wilcoxon rank sum test
wilcox.test(AUC[which(FLT3==0)],AUC[which(FLT3==1)], data=data)

lm.fit <- lm(AUC ~ FLT3, data=data)
summary(lm.fit)$coefficients[2,4]


#2. (i) Use a linear regression model to identify genes that are associated with AUC
# match the sample names
gex <- gex[,match(data$sample, colnames(gex))]
# Perform the linear regression
results <- apply(gex, 1, function(x) summary(lm(AUC ~ x, data=data))$coefficients[2,c(1,4)])
# Create a data frame with the results
results_df <- data.frame(Gene = rownames(gex),Estimate = results[1, ],p_value = results[2, ])
# Correct for multiple testing
results_df$FDR <- p.adjust(results_df$p_value, method = "BH")

significant_genes <- results_df[results_df$FDR < 0.05, ]
num_significant_genes <- nrow(significant_genes)
num_significant_genes

par(mfrow=c(1,1))
hist(results_df$p_value, main="P-value distribution", xlab="P-value", ylab="Frequency", breaks = 10)

# Create a histogram of p-values
hist(results_df$p_value, breaks = 50, main = "Distribution of p-values",
     xlab = "p-value", col = "lightblue", border = "black")
abline(h = length(results_df$p_value) / 50, col = "red", lty = 2)

# Select the gene with the strongest evidence of association
top_gene <- significant_genes[which.min(significant_genes$p_value),]

# Create the unadjusted model
unadjusted_model <- lm(AUC ~ gex[top_gene$Gene, ], data = data)
adjusted_model <- lm(AUC ~ gex[top_gene$Gene, ] + FLT3, data = data)

summary(unadjusted_model)
summary(adjusted_model)

#3. (i) cluster the samples based on gene expression
# Perform hierarchical clustering
normalized_gex <- scale(gex)

samplehc <- hclust(dist(t(normalized_gex)))
plot(samplehc)
# Cut the tree at level 5
clusters <- cutree(samplec, k = 5)

table(clusters)

# (ii) Draw a heatmap of the gene expression data
library(ComplexHeatmap)

# Dichotomize AUC into responder/non-responder
data$response <- ifelse(data$auc > median(data$auc), "Responder", "Non-responder")
#  Perform hierarchical clustering
samplehc <- hclust(dist(t(normalized_gex)))
genehc <- hclust(dist(normalized_gex))
# Define clusters (5 clusters for individuals)
clusters <- cutree(samplehc, k = 5)

annotation <- HeatmapAnnotation(
  Cluster = factor(clusters),
  Response = data$response,
  col = list(
    Cluster = structure(1:5, names = levels(factor(clusters))),  # Color for clusters
    Response = c("Responder" = "lightblue", "Non-responder" = "pink")  # Colors for AUC response
  )
)

Heatmap(
  normalized_gex,
  name = "Gene Expression",
  cluster_rows = genehc,
  cluster_columns = samplehc,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = annotation,
  column_title = "Individuals",
  row_title = "Genes",
  heatmap_legend_param = list(title = "Z-score"),
  column_split = 5
)

# 4. (i) Perform a multiple linear regression

# Select features and response variable
clinical_features <- data[, c("responseToInductionTx", "consensusAMLFusions",
                                       "FLT3-ITD", "NPM1", "RUNX1", "ASXL1",
                                       "TP53", "DNMT3A", "CEBPA", "NRAS",
                                       "wbcCount", "albumin", "creatinine",
                                       "hematocrit", "hemoglobin")]

# Fit linear regression model
model <- lm(data$auc ~ ., data = clinical_features)

summary(model)

# Predictions
predictions <- predict(model)
data$auc[data$sample == "BA2409"]
predictions[which(data$sample == "BA2409")]

# Plot in-sample predictions vs. observed AUC
plot(AUC, predictions,
     xlab = "Observed AUC", ylab = "Predicted AUC",
     main = "Predicted vs Observed AUC")
abline(0, 1, col = "red")

par(mfrow = c(1, 2))

# Residuals vs Fitted plot
plot(fitted(model), resid(model),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted",
     col = "darkgray", pch = 16)
abline(h = 0, col = "red", lwd = 2)

# Q-Q plot of residuals
qqnorm(resid(model), main = "Q-Q Plot of Residuals")
qqline(resid(model), col = "red", lwd = 2)


# 4. (ii) predict correlation coefficient
cor(data$auc, predictions)

# 10-Fold Cross-Validation
makeFolds <- function(k = 10, n, seedVal = 123456) {
  set.seed(seedVal)
  s = trunc(n / k)
  foldId = rep(c(1:k), s + 1)[1:n]
  sampleIDx = sample(seq(n))
  kFolds = split(sampleIDx, foldId)
  return(kFolds)
}

# Prepare data
n <- nrow(clinical_features)
folds <- makeFolds(10, n)
cv_predictions <- NULL
cv_observed <- NULL

# Perform 10-Fold Cross-Validation
for (i in 1:length(folds)) {
  fi <- folds[[i]]
  # Split data into training and test sets
  train_data <- clinical_features[-fi, ]
  train_response <- data$auc[-fi]
  test_data <- clinical_features[fi, ]
  test_response <- data$auc[fi]

  # Fit linear regression model on training data
  fold_model <- lm(train_response ~ ., data = train_data)

  # Predict on test set
  fold_predictions <- predict(fold_model, newdata = test_data)

  # Store predictions and observed values
  cv_predictions <- c(cv_predictions, fold_predictions)
  cv_observed <- c(cv_observed, test_response)
}

cor(cv_observed, cv_predictions)
