#' @title Predictions for a new dataset using an existing bartbma object
#' 
#' @description This function produces predictions for a new dataset using a previously obtained bartBMA object.
#' @param object A bartBMA object obtained using the barBMA function.
#' @param newdata Covariate matrix for new dataset.
#' @export 
#' @return A vector of predictions for the new dataset.

predict_bartBMA<-function(object,newdata){
  #preds<-get_BART_BMA_test_predictions(newdata,object$bic,object$sumoftrees,object$y_minmax)
  #orig_preds<-preds[[1]]
  #class(orig_preds)<-"predict.bartBMA"
  #orig_preds
  
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    ret<-preds_bbma_lin_alg_outsamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                    nrow(object$test_data),object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,
                                    object$test_data
    )
  }else{if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    ret<-preds_bbma_lin_alg_insamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                   object$a,object$sigma,0,object$nu,
                                   object$lambda#diff_inital_resids
    )
    
  }else{
    #if test data included in call to object
    ret<-preds_bbma_lin_alg_outsamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                    nrow(newdata), object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,
                                    newdata
    )
  }}
  
  
  class(ret)<-"predict.bartBMA"
  ret
  
}