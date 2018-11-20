varImpScores<-function(object,...){
  #object will be bartBMA object.
  imp_vars2=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees)
  res<-apply(imp_vars2[[4]],2,sum)
  class(res)<-"varImpScores.bartBMA"  
  res
}



