review_weights = function(cv_sl) {
  "
  - Taken from:
    https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html
  - input object must be an output from SuperLearner
  "
  meta_weights = coef(cv_sl)
  means = colMeans(meta_weights)
  sds = apply(meta_weights, MARGIN = 2,  FUN = sd)
  mins = apply(meta_weights, MARGIN = 2, FUN = min)
  maxs = apply(meta_weights, MARGIN = 2, FUN = max)
  # Combine the stats into a single matrix.
  sl_stats = cbind("mean(weight)" = means, "sd" = sds, "min" = mins, "max" = maxs)
  # Sort by decreasing mean weight.
  sl_stats[order(sl_stats[, 1], decreasing = TRUE), ]
}


roc.helper <- function(predictions, labels){
  # returns coordinates for roc curve
  stopifnot(length(predictions) == length(labels))
  
  n  <- length(predictions)
  df <- cbind(predictions, labels)
  df <- df[order(df[,1], decreasing = TRUE),]
  
  roc.coord <- matrix(data=0, nrow=n, ncol=2)
  
  for (i in 1:n){
    l <- as.vector(df[i,])[2]
    
    roc.coord[i,(l+1)] <- 1
  }
  
  roc.coord <- apply(roc.coord, 2, cumsum)
  stopifnot(roc.coord[n, 2] == sum(labels))
  roc.coord[,1] <- roc.coord[,1]/roc.coord[n,1]
  roc.coord[,2] <- roc.coord[,2]/roc.coord[n,2]
  
  roc.coord <- data.frame(roc.coord)
  colnames(roc.coord) <- c('FPR', 'TPR')
  
  return (roc.coord)
}

plot.roc <- function(predictions, labels){
  roc.coord <- roc.helper(predictions, labels)
  n <- length(roc.coord[,1])
  
  library(ggplot2)
  ggplot(data=roc.coord, aes(FPR, TPR)) +
    geom_path() + 
    labs(x='False Positive Rate (FPR)', 
         y='True Positive Rate (TPR)') +
    ggtitle(paste('ROC Curve\nArea under ROC =', round(sum(roc.coord$TPR)/n, digits = 3))) + 
    theme(plot.title = element_text(hjust = 0.5))
}


get_sens_spec <- function(predictions, labels){
  distance <- function(point){
    return (sqrt(sum((point - c(1, 0))^2)))
  }
  
  roc.coord <- roc.helper(predictions, labels)
  all_distances <- apply(roc.coord, 1, distance)
  best_coord <- all_distances == max(all_distances)
  
}