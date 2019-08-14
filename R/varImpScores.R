#' @title Variable importances as defined by Hernandez et al. (2018)
#' 
#' @description This measure defines the importance of a variable as the model-probability weighted sum of the number of splits on the variable of interest, divided by the sum over all variables of such weighted counts of splits.
#' @param object A bartBMA object obtained using the barBMA function.
#' @export 
#' @return A vector of variable importances. The variables are ordered in the same order that they occur in columns of the input covariate matrix used to obtain the input bartBMA object.
varImpScores<-function(object,...){
  #object will be bartBMA object.
  imp_vars2=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees)
  res<-apply(imp_vars2[[4]],2,sum)
  #create varImpPlot command
  vIP<-rep(NA,length(res))
  total_weighted_var_counts<-sum(res)
  #now get variable inclusion probabilities
  vIP<-res/total_weighted_var_counts
  class(vIP)<-"varImpScores.bartBMA"
  vIP
}