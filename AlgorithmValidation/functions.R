#' Calculate the offset of a lindear discriminant analysis
#' @param ldaoutput: linear discriminant analysis results returned by lda()
#' @param dt: data frame containing variables used in LDA analysis
#' @param vars: numeric indices of variables
#' @return offset
ldaInt <- function(ldaoutput, dt, vars){
  pred <- sum(dt[10, vars]*ldaoutput$scaling)
  obs <- predict(ldaoutput, dt)$x[10]
  return(obs - pred)
}

cv<- function (before, after){
  diff <- scale(after - before)^2
  cv <- vapply(1:nrow(diff), function(x){
    return(sum(diff[x,]))
  }, FUN.VALUE = numeric(1))

  return(cv)
}

rcv <- function(before, after, root = FALSE){
  bands <- colnames(before)
  diff <- scale(after - before)
  max <- data.frame(
    lapply(1:ncol(before), function(x){
      vapply(1:nrow(before), function(i){
        return (max(before[i,x], after[i,x]))
        }, FUN.VALUE = numeric(1))
      })
  )

  colnames(max) <- bands
  if(root){
    max2 <- sqrt(scale(max)^2)
  } else {
    max2 <- scale(max)^2
  }

  dm <- diff/max
  rcv <- vapply(1:nrow(dm), function(x){
    return(sum(dm[x,]))
  }, FUN.VALUE = numeric(1))

  return(rcv)
}

ndi <- function (band1, band2){
  ndi <- (band1 - band2)/(band1 + band2)
  ndi <- as.data.frame(ndi)
  colnames(ndi) <- 'ndi'
  return(ndi)
}

#' calculate the LDA coeffecients for a given habitat and disturbance combination
#'
#' @param hab string used in regex expresssion to select habitat type.
#' @param dist string used in regex expression to select disturbance type.
#' @return A named list containing the LDA coefficients and the roc curves
#'  generated from performing linear discriminatn analysis on the IW and MAD
#'  algorithm data subsets defined by \code{dist} and \code{hab}.
#' @examples
#' lda_analysis('Forest', 'Bare')
lda_analysis <- function(hab, dist){
  # filter dataframes to records matching habitat and disturbance type
  iw <- filter(IW, grepl(hab, habitat, ignore.case = TRUE),
               grepl(dist, disturbance, ignore.case = TRUE))
  mad <- filter(MAD, grepl(hab, habitat, ignore.case = TRUE),
                grepl(dist, disturbance, ignore.case = TRUE))

  # perform linear discriminant analysis
  MADlda <- lda(formula = change ~ chi+v1+v2+v3+v4+v5+v6,
                data = filter(mad, Train == 1),
                CV = TRUE)

  IWlda <- lda(formula = change ~ cv_z+nbr_z+ndsi_z+ndvi_z+ndwi_z+rcvmax_z,
               data = filter(iw, Train == 1),
               CV = TRUE)

  IWpred <- predict(IWlda, newdata = filter(iw, Train == 0))
  MADpred <- predict(MADlda, newdata = filter(mad, Train == 0))

  IWroc <- roc(iw$change[iw$Train == 0], predictor = IWpred$x[,1])
  MADroc <- roc(mad$change[mad$Train == 0], predictor = MADpred$x[,1])
  return(list('IW'=IWroc, 'MAD'=MADroc, "IWcoeffs" = IWlda$scaling, "MADcoeffs" = MADlda$scaling))
}

#' calculate the QDA coeffecients for a given habitat and disturbance combination
#'
#' @param hab string used in regex expresssion to select habitat type.
#' @param dist string used in regex expression to select disturbance type.
#' @return A named list containing the LDA coefficients and the roc curves
#'  generated from performing linear discriminatn analysis on the IW and MAD
#'  algorithm data subsets defined by \code{dist} and \code{hab}.
#' @examples
#' lda_analysis('Forest', 'Bare')
qda_analysis <- function(hab, dist){
  # filter dataframes to records matching habitat and disturbance type
  iw <- filter(IW, grepl(hab, habitat, ignore.case = TRUE),
               grepl(dist, disturbance, ignore.case = TRUE))
  mad <- filter(MAD, grepl(hab, habitat, ignore.case = TRUE),
                grepl(dist, disturbance, ignore.case = TRUE))

  # perform quadraic discriminant analysis
  MADqda <- qda(formula = change ~ chi+v1+v2+v3+v4+v5+v6,
                data = filter(mad, Train == 1))

  IWqda <- qda(formula = change ~ cv_z+nbr_z+ndsi_z+ndvi_z+ndwi_z+rcvmax_z,
               data = filter(iw, Train == 1))

  IWpred <- predict(IWqda, newdata = filter(iw, Train == 0))
  MADpred <- predict(MADqda, newdata = filter(mad, Train == 0))

  IWroc <- roc(iw$change[iw$Train == 0], predictor = IWpred$x[,1])
  MADroc <- roc(mad$change[mad$Train == 0], predictor = MADpred$x[,1])
  return(list('IW'=IWroc, 'MAD'=MADroc, "IWcoeffs" = IWlda$scaling, "MADcoeffs" = MADlda$scaling))
}

#' Generate plotly graph displaying ROC curve for LDA output of IW and MAD data
#'
#' Designed to be used in conjunction with \code{lda_analysis} function, taking the
#' \code{IW} and \code{MAD} named output.
#'
#' @param iwroc object returned from \code{roc} function using IW data as input
#' @param madroc object returned form \code{roc} function using MAD data as input
#' @return plotly graph object with two line traces
#' @examples
#' ldaOutput <- lda_analysis("Forest", 'Bare')
#' plot_ROC_crv(ladOutput$IW, ldaOutput$MAD)
#'
plot_ROC_crv <- function(iwroc, madroc){
  plot_ly(type = 'scatter', mode = 'lines')%>%
    add_trace(name = "MAD", x = ~1-madroc$specificities,
              y = ~madroc$sensitivities, line = list(color = 'blue'),
              hoverinfo = 'text',
              text = ~paste('LDA score:', madroc$thresholds,
                            "<br>True Positive:", madroc$sensitivities,
                            "<br>False Positive:", 1-madroc$specificities,
                            sep = " "))%>%
    add_trace(name = "LCC", x = ~1-iwroc$specificities,
              y = ~iwroc$sensitivities, line = list(color = 'orange'),
              hoverinfo = 'text',
              text = ~paste('LDA score:', iwroc$thresholds,
                            "<br>True Positive:", iwroc$sensitivities,
                            "<br>False Positive:", 1-iwroc$specificities,
                            sep = " "))%>%
    layout(title = "ROC: All Changes",
           legend = list(x = 0.75, y = 0.2),
           xaxis = list(title = "False Positive Rate"),
           yaxis = list(title = "True Positive Rate"))
}
