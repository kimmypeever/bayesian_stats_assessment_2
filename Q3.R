# ------------------------------
# Q3 Applying Bayesian Networks
# ------------------------------

# Load required packages
library(bnlearn)
library(dplyr)
library(Rgraphviz)

# Load the dataset
nursery <- read.csv("nursery.data", header = FALSE, stringsAsFactors = TRUE)

# Assign column names (from UCI dataset description)
colnames(nursery) <- c(
  "parents", "has_nurs", "form", "children", "housing",
  "finance", "social", "health", "class"
)

# View structure of the dataset
str(nursery)

# Check for missing values
sum(is.na(nursery))

# Confirm all variables are factors (categorical)
sapply(nursery, is.factor)

# Display number of levels per variable
sapply(nursery, nlevels)

# Summary of the target variable (class imbalance check)
table(nursery$class)
prop.table(table(nursery$class))

# Basic descriptive overview
summary(nursery)

# Visualize class distribution
barplot(table(nursery$class),
        main = "Distribution of Target Classes",
        xlab = "Class",
        ylab = "Count",
        col = "lightblue")
#-----------------------------------------------------
# Learn the structure using the Hill-Climbing algorithm
bn_model <- hc(nursery)

# Plot with Graphviz
graphviz.plot(bn_model, main = "Bayesian Network Structure for Nursery Data")

# View the structure (edges)
bn_model

# Fit the parameters (conditional probability tables)
fitted_bn <- bn.fit(bn_model, nursery)

# Inspect the CPTs for the class variable
fitted_bn$class


#---------------------------------------------------------
# Perform 10-fold cross-validation using log-likelihood as the scoring metric
set.seed(8)
cv_results <- bn.cv(nursery, bn_model, loss = "logl", k = 10)

# View cross-validation results
cv_results

#---------------------------------------------------------
# Priority recommendation given slightly problematic social 
# conditions and a recommended health status
cpquery(
  fitted = fitted_bn,
  event = (class == "priority"),
  evidence = list(social = "slightly_prob", health = "recommended"),
  method = "lw",
  n = 10000
)


# Very recommended outcome for applicants with non-problematic social
# conditions and convenient financial status
cpquery(
  fitted = fitted_bn,
  event = (class == "very_recom"),
  evidence = list(social = "nonprob", finance = "convenient"),
  method = "lw",
  n = 10000
)


# Not recommended outcome for applicants with problematic
# social situations and poor health
cpquery(
  fitted = fitted_bn,
  event = (class == "not_recom"),
  evidence = list(social = "problematic", health = "not_recom"),
  method = "lw",
  n = 10000
)

# Overall class distribution for applicants with convenient housing conditions
table(cpdist(
  fitted_bn,
  nodes = "class",
  evidence = (housing == "convenient")
)) / 10000


