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

plot.roc.conf.int <- function(predictions, labels, bs_predictions){
  roc.coord <- roc.helper(predictions, labels)
  n <- length(roc.coord[,1])
  
  fpr <- roc.coord$FPR
  tpr <- c()
  q <- c(0, 0)
  auroc.vals <- c()
  
  for (i in 1:dim(bs_predictions)[2]){
    temp.coord <- roc.helper(bs_predictions[,i], labels)
    
    temp_tpr <- c()
    
    for (f in fpr){
      temp_tpr <- c(temp_tpr, max(roc.coord[(abs(temp.coord[,1] - f) < 0.0001),]$TPR))
    }
    
    temp_tpr[temp_tpr == -Inf] = 0
    tpr <- cbind(tpr, temp.coord$TPR)
    auroc.vals <- c(sum(temp.coord$TPR)/n, auroc.vals)
  }
  
  auroc.confint <- round(as.vector(quantile(auroc.vals, c(0.025, 0.975))), digit=2)
  
  for (i in 1:n){
    q <- rbind(q, quantile(tpr[i,], probs=c(0.025, 0.975))) # 95% conf int
  }
  
  q <- q[-c(1),]
  roc.coord$lower <- q[,1]
  roc.coord$upper <- q[,2]
  
  library(ggplot2)
  ggplot(data=roc.coord, aes(FPR, TPR)) +
    geom_path(aes(y=TPR)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='red', alpha='0.5') +
    labs(x='False Positive Rate (FPR)', 
         y='True Positive Rate (TPR)') +
    ggtitle(paste0('ROC Curve\nArea under ROC = ', round(sum(roc.coord$TPR)/n, digits = 2), ' (', auroc.confint[1], ', ', auroc.confint[2], ')')) + 
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