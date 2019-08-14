#' @title Predictions for a new dataset using an existing probit_bartBMA object
#' 
#' @description This function produces predictions for a new dataset using a previously obtained bartBMA object.
#' @param object A probit_bartBMA object obtained using the probit_bartBMA function.
#' @param newdata Covariate matrix for new dataset.
#' @export 
#' @return A vector of predictions for the new dataset.

predict_probit_bartBMA<-function(object,newdata){
  
  
  preds<-get_BART_BMA_test_predictions(newdata,object$bic,object$sumoftrees,object$y_minmax)
  
  ret <- list()
  
  probs <-pnorm(preds[[1]])
  pred_binary <- ifelse(probs<=0.5,0,1)
  
  ret$probs <- probs
  ret$pred_binary <- pred_binary
  
  ret
}