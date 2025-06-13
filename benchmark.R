data<-read.table("PNF_EAD.csv",sep=",",h=T,na.string="NA")
library(dplyr)


data%>%select(outcome_PNF_EAD,BILI,PLQ,TEMPS_INTERVENTION,REPERFUSION_STEATOSE_MICRO)->small
small%>%dplyr::rename(REP_MICROSTEAT="REPERFUSION_STEATOSE_MICRO")->small

##install.packages("naniar")  library(naniar)

# Heatmap missing data
gg_miss_var(small)
gg_miss_upset(small)
vis_miss(small)






library(naniar)

# Test MCAR on dataframe 
mcar_test(small)






#############
## imputation missforest

library(missForest)
df_impute <- missForest(small)

# print results
print(df_impute$ximp)  # imputed data
print(df_impute$OOBerror)  # error


save(df_impute,file="imputed.rda")
load("imputed.rda")
library(rms)

m <- lrm(outcome_PNF_EAD ~ BILI + PLQ +TEMPS_INTERVENTION + REPERFUSION_STEATOSE_MICRO, 
             data = df, x = TRUE, y = TRUE)
m
######
df<-as.data.frame(df_impute$ximp)

df$outcome_PNF_EAD<-as.factor(df$outcome_PNF_EAD)

df%>%dplyr::rename(Class="outcome_PNF_EAD")->df

data<-df
data$Class <- factor(data$Class, levels = c(0,1), labels = c("Class0", "Class1"))
colnames(data) <- make.names(colnames(data))

# load packages
library(caret)
library(glmnet)


# split cohort
set.seed(45678)
trainIndex <- createDataPartition(data$Class, p = 0.6, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# list of models
models <- list(
  "glm" = "glm",              # Régression logistique
  "rpart" = "rpart",          # Arbre de décision
  "knn" = "knn",              # K-plus proches voisins
  "lasso" = "glmnet",         # Régression Lasso
  "ridge" = "glmnet",         # Régression Ridge
  "elasticnet" = "glmnet"     # Elastic Net
)

# save results
results <- list()

# loop benchmark models
for (model_name in names(models)) {
  cat("Entraînement du modèle :", model_name, "\n")
  
  control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  
  if (model_name %in% c("lasso", "ridge", "elasticnet")) {
    #  glmnet parameters
    fit <- train(
      Class ~ ., data = trainData, method = "glmnet",
      tuneGrid = expand.grid(alpha = c(1, 0, 0.5), lambda = seq(0.01, 0.1, 0.01)),
      metric = "ROC", trControl = control
    )

  } else {
    # other models
    fit <- train(Class ~ ., data = trainData, method = models[[model_name]], metric = "ROC", trControl = control)
  }
  
  # accuracy
  predictions <- predict(fit, newdata = testData)
  accuracy <- sum(predictions == testData$Class) / nrow(testData)
  
  # save
  results[[model_name]] <- list(Model = fit, Accuracy = accuracy)
}

# print results
for (model_name in names(results)) {
  cat("Modèle :", model_name, "- Exactitude :", results[[model_name]]$Accuracy, "\n")
}


# load ggplot2
library(ggplot2)

# dataframe results
results_df <- data.frame(
  Model = character(),
  Accuracy = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  stringsAsFactors = FALSE
)

# CI95
for (model_name in names(results)) {
  accuracy <- results[[model_name]]$Accuracy
  n <- nrow(testData)
  se <- sqrt((accuracy * (1 - accuracy)) / n) 
  z <- 1.96  # Pour un intervalle de confiance à 95 %
  lower_ci <- accuracy - z * se
  upper_ci <- accuracy + z * se
  
  results_df <- rbind(
    results_df,
    data.frame(
      Model = model_name,
      Accuracy = accuracy,
      Lower_CI = max(0, lower_ci), 
      Upper_CI = min(1, upper_ci)   
    )
  )
}



# order
results_df <- results_df[order(results_df$Accuracy, decreasing = TRUE), ]

# Forest Plot with ggplot2
ggplot(results_df, aes(x = Accuracy, y = reorder(Model, Accuracy))) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.3, color = "darkgray") +
  geom_text(aes(label = round(Accuracy, 3)), hjust = -0.5,vjust=-0.5, color = "black") +
  theme_minimal() +
  labs(
    title = "Forest Plot of models with CI95 - cohort split 0.6/0.4",
    x = "(Accuracy)",
    y = "Models"
  ) +
  theme(plot.title = element_text(hjust = 0.5))+coord_cartesian(xlim = c(0.85, 0.95))









