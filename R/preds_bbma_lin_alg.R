#' @title Predictions for bart-bma output obtained from the posterior probability weighted averaged of the posterior means for each model
#' 
#' @description This function produces predictions from BART-BMA by obtaining the posterior probability weighted averaged of the posterior means for each model.
#' @param object bartBMA object obtained from function bartBMA
#' @param num_iter Total number of iterations of the Gibbs sampler (including burn-in).
#' @param burnin Number of burn-on iterations of the Gibbs sampler.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residua;s.
#' @export 
#' @return The output is a list of length one. The one element in this list is a vector of prediction intervals???

preds_bbma_lin_alg<-function(object,num_iter,burnin,newdata=NULL,update_resids=1,trainingdata){
  #object will be bartBMA object.
  

  
  
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

  

  ret
  
}