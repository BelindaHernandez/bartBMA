predict.bartBMA<-function(object,newdata){
  preds<-get_BART_BMA_test_predictions(newdata,object$bic,object$sumoftrees,object$y_minmax)
  orig_preds<-preds[[1]]
  class(orig_preds)<-"predict.bartBMA"
  orig_preds
}