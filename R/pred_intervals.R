pred_intervals<-function(object,num_iter,burnin,l_quant,u_quant,newdata=NULL){
  if(l_quant>0.5 ||u_quant<0 ||u_quant>1){stop("Lower quantile must be lower than 0.5 and greater than 0")}
  if(u_quant<0.5 ||u_quant<0 ||u_quant>1){stop("Upper quantile must be greater than 0.5 and less than 1")}
  #object will be bartBMA object.
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    gs_chains<-gibbs_sampler(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                             nrow(object$test_data),object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals,object$test_data)
  }else if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    gs_chains<-gibbs_sampler2(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                              object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals)
    
  }else{
    #if test data included in call to object
    gs_chains<-gibbs_sampler(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,num_iter, burnin,object$nrowTrain,
                             nrow(newdata), object$a,object$sigma,0,object$nu,object$lambda,object$sum_residuals,newdata)
  }
  
  y_posterior_sum_trees<-gs_chains[[4]]
  y_orig_post_sum_trees<-gs_chains[[5]]
  sum_of_tree_BIC<-unlist(object$bic)
  weights<-exp(sum_of_tree_BIC-(max(sum_of_tree_BIC)+log(sum(exp(sum_of_tree_BIC-max(sum_of_tree_BIC))))))
  sigma_chains<-gs_chains[[7]]
  final_length=num_iter-burnin
  num_its_to_sample<-round(weights*final_length)
  final_sigma_chain<-numeric(0)
  
  final_y_chain<-matrix(nrow=0,ncol=ncol(y_posterior_sum_trees[[1]]))
  final_yorig_chain<-matrix(nrow=0,ncol=ncol(y_orig_post_sum_trees[[1]]))
  
  for(i in 1:length(sigma_chains)){
    sample_its<-sample(burnin:num_iter,num_its_to_sample[i])
    final_sigma_chain<-c(final_sigma_chain,sigma_chains[[i]][sample_its])
    #now do the same for predicted response updates
    post_y_i<-y_posterior_sum_trees[[i]]
    post_yorig_i<-y_orig_post_sum_trees[[i]]
    final_y_chain<-rbind(final_y_chain,post_y_i[sample_its,])
    final_yorig_chain<-rbind(final_yorig_chain,post_yorig_i[sample_its,])
    PI<-apply(final_yorig_chain,2,function(x)quantile(x,probs=c(l_quant,0.5,u_quant)))
  }
  
  
  ret<-list()
  length(ret)<-1
  ret[[1]]<-PI
  
  class(ret)<-"pred_intervals.bartBMA"  
  ret
}