library(ggplot2)
library(ranger)
library(fastshap)

# Read in the data and clean it up a bit
titanic <- titanic::titanic_train
features <- c(
  "Survived",  # passenger survival indicator
  "Pclass",    # passenger class
  "Sex",       # gender
  "Age",       # age
  "SibSp",     # number of siblings/spouses aboard
  "Parch",     # number of parents/children aboard
  "Fare",      # passenger fare
  "Embarked"   # port of embarkation
)
titanic <- titanic[, features]
titanic$Survived <- as.factor(titanic$Survived)
titanic <- na.omit(titanic)  # ...umm?

set.seed(1046)  # for reproducibility
rfo <- ranger(Survived ~ ., data = titanic, probability = TRUE)
# Prediction wrapper for `fastshap::explain()`; has to return a 
# single (atomic) vector of predictions
pfun <- function(object, newdata) {  # computes prob(Survived=1|x)
  predict(object, data = newdata)$predictions[, 2]
}

# Estimate feature contributions for each imputed training set
X <- subset(titanic, select = -Survived)  # features only!

set.seed(1051)  # for reproducibility
(ex.all <- explain(rfo, X = X, nsim = 100, adjust = TRUE,  pred_wrapper = pfun))

vi_shap(rfo, train = X, pred_wrapper = pfun, nsim = 10)

vi_shap(stack_list[[1]], train = feats_list %>% bind_rows %>% select(-gid, -subtype, -label, -src), feature_names = c("A_Bias_HA"),  pred_wrapper = pfun, nsim = 10)
