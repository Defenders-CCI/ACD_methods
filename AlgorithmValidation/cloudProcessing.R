library(caret)
library(dplyr)
library(MASS)
library(mvtnorm)
library(plotly)
library(pROC)
library(stringr)
library(tidyr)

cloud <- data.frame()

for (file in list.files("ChangeValidation/data/clouds", pattern = '.csv')){
  dat <- read.csv(paste("ChangeValidation/data/clouds", file, sep = "/"), header = TRUE)
  colnames(dat) <- trimws(
    tolower(
      gsub("\\.", "", colnames(dat))
    )
  )
  dat$Class <- paste(str_split_fixed(dat$systemindex, "_", n = 3)[,1], str_split_fixed(dat$systemindex, "_", n = 3)[,2], sep = "")
  dat$geo <- NULL
  dat$Loc <- strsplit(file, "_")[[1]][1]
  cloud <- bind_rows(cloud, dat)
}

cloud$Class <- vapply(cloud$Class, function(x){
  switch(x, '11' = 'clear', '12' = 'cloud', '20' = 'shadow')},
  FUN.VALUE = character(1), USE.NAMES = FALSE)

cloud$train = rbinom(nrow(cloud), 1, 0.7)
cloud$shadow <- ifelse(cloud$Class == "shadow", 1, 0)
cloud$cloud <- ifelse(cloud$Class == 'cloud', 1, 0)

shadow_lda <- lda(formula = shadow ~ b11+b12+b8+c1+c2+c3, data = filter(cloud, train == 1))
shadow_pred <- predict(shadow_lda, newdata = filter(cloud, train == 0))
shadow_roc <- roc(response = cloud$shadow[cloud$train == 0], predictor = shadow_pred$x[,1])

shadowtree <- rpart(shadow ~ b11+b12+b8+c1+c2+c3, data = filter(cloud, train == 1), method = 'class')
shadowtree_pred <- predict(shadowtree, newdata = filter(cloud, train == 0))
shadowtree_roc <- roc(cloud$shadow[cloud$train == 0], predictor = shadowtree_pred[,2])

cloud_lda <- lda(formula = cloud ~ ndmi_z + ndsi_z + vis_z + cirrus_z + cloudscore_z, data = cloud)
cloud_pred <- predict(cloud_lda, newdata = cloud)
cloud_roc <- roc(response = cloud$cloud, predictor = cloud_pred$x[,1])

cloudtree <- rpart(cloud ~ ndmi_z + ndsi_z + vis_z + cirrus_z + cloudscore_z, data = cloud, method = 'class')
cloudtree_pred <- predict(cloudtree, newdata = cloud)
cloudtree_roc <- roc(cloud$cloud, predictor = cloudtree_pred[,2])


plot_ly(type = 'scatter', mode = 'lines')%>%
  add_trace(name = "Clouds Tree", x = ~1-cloudtree_roc$specificities,
            y = ~cloudtree_roc$sensitivities, line = list(color = 'purple'),
            hoverinfo = 'text',
            text = ~paste('LDA Score:', cloudtree_roc$thresholds,
                          "<br>True Positive:", cloudtree_roc$sensitivities,
                          "<br>False Positive:", 1-cloudtree_roc$specificities,
                          sep = " "))%>%
  add_trace(name = "Clouds LDA", x = ~1-cloud_roc$specificities,
            y = ~cloud_roc$sensitivities, line = list(color = 'purple', dash = 'dash'),
            hoverinfo = 'text',
            text = ~paste('LDA Score:', cloud_roc$thresholds,
                          "<br>True Positive:", cloud_roc$sensitivities,
                          "<br>False Positive:", 1-cloud_roc$specificities,
                          sep = " "))%>%
  add_trace(name = "Shadow LDA", x = ~1-shadow_roc$specificities,
            y = ~shadow_roc$sensitivities, line = list(color = 'green'),
            hoverinfo = 'text',
            text = ~paste('LDA Score:', shadow_roc$thresholds,
                          '<br>True Positive:', shadow_roc$sensitivities,
                          '<br>False Positive:',
                          1-shadow_roc$specificities))%>%
  add_trace(name = "Shadow Tree", x = ~1-shadowtree_roc$specificities,
            y = ~shadowtree_roc$sensitivities, line = list(color = 'green', dash = 'dash'),
            hoverinfo = 'text',
            text = ~paste('LDA Score:', shadowtree_roc$thresholds,
                          '<br>True Positive:', shadowtree_roc$sensitivities,
                          '<br>False Positive:',
                          1-shadowtree_roc$specificities))%>%
  layout(title = "ROC: All Changes",
         legend = list(x = 0.75, y = 0.2),
         xaxis = list(title = "False Positive Rate"),
         yaxis = list(title = "True Positive Rate"))

cloud$train <- rbinom(nrow(cloud), 1, 0.8)
