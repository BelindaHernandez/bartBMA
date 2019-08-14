#' @title Predictions for a new dataset using an existing bartbma object
#' 
#' @description This function produces predictions for a new dataset using a previously obtained bartBMA object.
#' @param object A bartBMA object obtained using the barBMA function.
#' @param newdata Covariate matrix for new dataset.
#' @export 
#' @return A vector of predictions for the new dataset.

predict_bartBMA<-function(object,newdata){
  preds<-get_BART_BMA_test_predictions(newdata,object$bic,object$sumoftrees,object$y_minmax)
  orig_preds<-preds[[1]]
  class(orig_preds)<-"predict.bartBMA"
  orig_preds
}