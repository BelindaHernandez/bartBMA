#' @title Predictions for a new dataset using an existing probit_bartBMA object
#' 
#' @description This function produces predictions for a new dataset using a previously obtained bartBMA object.
#' @param object A probit_bartBMA object obtained using the probit_bartBMA function.
#' @param newdata Covariate matrix for new dataset.
#' @export 
#' @return A vector of predictions for the new dataset.

predict_probit_bartBMA<-function(object,newdata){
  
  
  #preds<-get_BART_BMA_test_predictions(newdata,object$bic,object$sumoftrees,object$y_minmax)
  
  
  
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    preds<-preds_bbma_lin_alg_outsamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                    nrow(object$test_data),object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,
                                    object$test_data
    )
  }else{if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    preds<-preds_bbma_lin_alg_insamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                   object$a,object$sigma,0,object$nu,
                                   object$lambda#diff_inital_resids
    )
    
  }else{
    #if test data included in call to object
    preds<-preds_bbma_lin_alg_outsamp(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                                    nrow(newdata), object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,
                                    newdata
    )
  }}
  
  
  
  
  
  
  
  
  ret <- list()
  
  #probs <-pnorm(preds[[1]])
  probs <-pnorm(preds)
  
  pred_binary <- ifelse(probs<=0.5,0,1)
  
  ret$probs <- probs
  ret$pred_binary <- pred_binary
  
  ret
}