#' @title Obtains ITE Predictions (in-sample)  and Prediction intervals using bartBMA, a Gibbs sampler, and the method described by Hill (2011)
#' 
#' @description This function produces predictions for a new dataset using a previously obtained bartBMA object.
#' @param object bartBMA object obtained from function bartBMA
#' @param num_iter Total number of iterations of the Gibbs sampler (including burn-in).
#' @param burnin Number of burn-on iterations of the Gibbs sampler.
#' @param l_quant Lower quartile of the prediction interval.
#' @param u_quant Upper quartile of the prediction interval.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residua;s.
#' @param x_covariates Covaraite matrix for training bartBMA.
#' @param z_train treatment vector for traiing bartBMA.
#' @param y_train outcome vector for training bartBMA.
#' @param a This is a parameter that influences the variance of terminal node parameter values. Default value a=3.
#' @param nu This is a hyperparameter in the distribution of the variance of the error term. THe inverse of the variance is distributed as Gamma (nu/2, nu*lambda/2). Default value nu=3.
#' @param sigquant ??
#' @param c This determines the size of Occam's Window
#' @param pen This is a parameter used by the Pruned Exact Linear Time Algorithm when finding changepoints. Default value pen=12.
#' @param num_cp This is a number between 0 and 100 that determines the proportion of changepoints proposed by the changepoint detection algorithm to keep when growing trees. Default num_cp=20.
#' @param x.test Test data covariate matrix. Default x.test=matrix(0.0,0,0).
#' @param num_rounds Number of trees. (Maximum number of trees in a sum-of-tree model). Default num_rounds=5.
#' @param alpha Parameter in prior probability of tree node splitting. Default alpha=0.95
#' @param beta Parameter in prior probability of tree node splitting. Default beta=1
#' @param split_rule_node Binary variable. If equals 1, then find a new set of potential splitting points via a changepoint algorithm after adding each split to a tree. If equals zero, use the same set of potential split points for all splits in a tree. Default split_rule_node=0.
#' @param gridpoint Binary variable. If equals 1, then a grid search changepoint detection algorithm will be used. If equals 0, then the Pruned Exact Linear Time (PELT) changepoint detection algorithm will be used (Killick et al. 2012). Default gridpoint=0.
#' @param maxOWsize Maximum number of models to keep in Occam's window. Default maxOWsize=100.
#' @param num_splits Maximum number of splits in a tree
#' @param gridsize This integer determines the size of the grid across which to search if gridpoint=1 when finding changepoints for constructing trees.
#' @param zero_split Binary variable. If equals 1, then zero split trees can be included in a sum-of-trees model. If equals zero, then only trees with at least one split can be included in a sum-of-trees model.
#' @param only_max_num_trees Binary variable. If equals 1, then only sum-of-trees models containing the maximum number of trees, num_rounds, are selected. If equals 0, then sum-of-trees models containing less than num_rounds trees can be selected. The default is only_max_num_trees=1.
#' @param min_num_obs_for_split This integer determines the minimum number of observations in a (parent) tree node for the algorithm to consider potential splits of the node.
#' @param min_num_obs_after_split This integer determines the minimum number of observations in a child node resulting from a split in order for a split to occur. If the left or right chikd node has less than this number of observations, then the split can not occur.
#' @param exact_residuals Binary variable. If equal to 1, then trees are added to sum-of-tree models within each round of the algorithm by detecting changepoints in the exact residuals. If equals zero, then changepoints are detected in residuals that are constructed from approximate predictions.
#' @param spike_tree If equal to 1, then the Spike-and-Tree prior will be used, otherwise the standard BART prior will be used. The number of splitting variables has a beta-binomial prior. The number of terminal nodes has a truncated Poisson prior, and then a uniform prior is placed on the set of valid constructions of trees given the splitting variables and number of terminal nodes.
#' @param lambda_poisson This is a parameter for the Spike-and-Tree prior. It is the parameter for the (truncated and conditional on the number of splitting variables) Poisson prior on the number of terminal nodes.
#' @export 
#' @return A vector of Individual Treatment Effect Estimates.

GS_ITEs_PIs_bartBMA_new_initials<-function(num_iter,burnin,l_quant,u_quant,newdata=NULL,update_resids=1,
                              x_covariates,z_train ,y_train,
                              a1=3,nu1=3,sigquant1=0.9,c1=1000,
                              pen1=12,num_cp1=20,x.test1=matrix(0.0,0,0),
                              num_rounds1=5,alpha1=0.95,beta1=2,split_rule_node1=0,
                              gridpoint1=0,maxOWsize1=100,num_splits1=5,gridsize1=10,zero_split1=1,only_max_num_trees1=1,
                              min_num_obs_for_split1=2, min_num_obs_after_split1=2,
                              exact_residuals1=1,spike_tree1=0,lambda_poisson1=10, less_greedy1=0){
  
  
  
 
  
  
  
  x_train <- cbind(z_train,x_covariates)
  
  trained_bart_BMA <- bartBMA(x.train= x_train,y.train =y_train ,
                              a=a1,nu=nu1,sigquant=sigquant1,c=c1,
                              pen=pen1,num_cp=num_cp1,x.test=x.test1,
                              num_rounds=num_rounds1,alpha=alpha1,beta=beta1,split_rule_node=split_rule_node1,
                              gridpoint=gridpoint1,maxOWsize=maxOWsize1,num_splits=num_splits1,gridsize=gridsize1,zero_split=zero_split1,only_max_num_trees=only_max_num_trees1,
                              min_num_obs_for_split=min_num_obs_for_split1, min_num_obs_after_split=min_num_obs_after_split1,
                              exact_residuals=exact_residuals1,spike_tree=spike_tree1,lambda_poisson=lambda_poisson1,less_greedy1=less_greedy)
  
  
  scaled_train_y <- scale_response(min(trained_bart_BMA$response),max(trained_bart_BMA$response),-0.5,0.5,trained_bart_BMA$response)
  
  # diff_inital_resids <- list()
  # resid_length <- length(trained_bart_BMA$sum_residuals[[1]][[1]])
  # #initial_partial_resids <- rep(0.8*mean(scaled_train_y),resid_length )
  # 
  # 
  # for(i in 1:length(trained_bart_BMA$sum_residuals)){
  #   diff_inital_resids[[i]] <- rep(list(scaled_train_y- ((length(trained_bart_BMA$sum_residuals[[i]])-1 )/ length(trained_bart_BMA$sum_residuals[[i]]))*mean(scaled_train_y)), 
  #                                  length(trained_bart_BMA$sum_residuals[[i]])) 
  # }
  
  get_resids <- get_initial_resids(x_train,trained_bart_BMA$sumoftrees,scaled_train_y)
  diff_inital_resids<- get_resids[[1]]
  new_pred_list1 <- get_resids[[2]]
  
  
  if(update_resids==0){
    if(is.null(newdata) && length(trained_bart_BMA)==16){
      all_treated_data <- cbind(rep(1,nrow(trained_bart_BMA$test_data)), trained_bart_BMA$test_data)
      all_control_data <- cbind(rep(0,nrow(trained_bart_BMA$test_data)), trained_bart_BMA$test_data)      
      #if test data specified separately
      gs_chains<-gibbs_sampler_ITE_no_update_new_inits(trained_bart_BMA$sumoftrees,trained_bart_BMA$obs_to_termNodesMatrix,trained_bart_BMA$response,trained_bart_BMA$bic,num_iter, burnin,trained_bart_BMA$nrowTrain,
                                             nrow(trained_bart_BMA$test_data),trained_bart_BMA$a,trained_bart_BMA$sigma,0,
                                             trained_bart_BMA$nu,trained_bart_BMA$lambda,
                                             diff_inital_resids,
                                             all_treated_data,all_control_data,
                                             new_pred_list1)
    }else if(is.null(newdata) && length(trained_bart_BMA)==14){
      all_treated_data <- cbind(rep(1,nrow(x_covariates)), x_covariates)
      all_control_data <- cbind(rep(0,nrow(x_covariates)), x_covariates)
      #else return Pred Ints for training data
      gs_chains<-gibbs_sampler_ITE_no_update2_new_inits(trained_bart_BMA$sumoftrees,trained_bart_BMA$obs_to_termNodesMatrix,trained_bart_BMA$response,trained_bart_BMA$bic,num_iter, burnin,trained_bart_BMA$nrowTrain,
                                              trained_bart_BMA$a,trained_bart_BMA$sigma,0,trained_bart_BMA$nu,trained_bart_BMA$lambda,
                                              diff_inital_resids,
                                              all_treated_data,all_control_data,
                                              new_pred_list1)
      
    }else{
      
      all_treated_data <- cbind(rep(1,nrow(newdata)), newdata)
      all_control_data <- cbind(rep(0,nrow(newdata)), newdata)
      #if test data included in call to object
      gs_chains<-gibbs_sampler_ITE_no_update_new_inits(trained_bart_BMA$sumoftrees,trained_bart_BMA$obs_to_termNodesMatrix,trained_bart_BMA$response,trained_bart_BMA$bic,num_iter, burnin,trained_bart_BMA$nrowTrain,
                                             nrow(newdata), trained_bart_BMA$a,trained_bart_BMA$sigma,0,trained_bart_BMA$nu,trained_bart_BMA$lambda,
                                             diff_inital_resids,
                                             all_treated_data,all_control_data,
                                             new_pred_list1)
    }
  }else{
    if(is.null(newdata) && length(trained_bart_BMA)==16){
      all_treated_data <- cbind(rep(1,nrow(trained_bart_BMA$test_data)), trained_bart_BMA$test_data)
      all_control_data <- cbind(rep(0,nrow(trained_bart_BMA$test_data)), trained_bart_BMA$test_data)      
      
      #if test data specified separately
      gs_chains<-gibbs_sampler_ITE_new_inits(trained_bart_BMA$sumoftrees,trained_bart_BMA$obs_to_termNodesMatrix,trained_bart_BMA$response,trained_bart_BMA$bic,num_iter, burnin,trained_bart_BMA$nrowTrain,
                                   nrow(trained_bart_BMA$test_data),trained_bart_BMA$a,trained_bart_BMA$sigma,0,
                                   trained_bart_BMA$nu,trained_bart_BMA$lambda,
                                   diff_inital_resids,
                                   all_treated_data,all_control_data,
                                   new_pred_list1)
      
    }else if(is.null(newdata) && length(trained_bart_BMA)==14){
      all_treated_data <- cbind(rep(1,nrow(x_covariates)), x_covariates)
      all_control_data <- cbind(rep(0,nrow(x_covariates)), x_covariates)
      #else return Pred Ints for training data
      gs_chains<-gibbs_sampler_ITE2_new_inits(trained_bart_BMA$sumoftrees,trained_bart_BMA$obs_to_termNodesMatrix,trained_bart_BMA$response,trained_bart_BMA$bic,num_iter, burnin,trained_bart_BMA$nrowTrain,
                                    trained_bart_BMA$a,trained_bart_BMA$sigma,0,trained_bart_BMA$nu,trained_bart_BMA$lambda,
                                    diff_inital_resids,
                                    all_treated_data,all_control_data,
                                    new_pred_list1)
      
    }else{
      
      all_treated_data <- cbind(rep(1,nrow(newdata)), newdata)
      all_control_data <- cbind(rep(0,nrow(newdata)), newdata)
      #if test data included in call to object
      gs_chains<-gibbs_sampler_ITE_new_inits(trained_bart_BMA$sumoftrees,trained_bart_BMA$obs_to_termNodesMatrix,trained_bart_BMA$response,trained_bart_BMA$bic,num_iter, burnin,trained_bart_BMA$nrowTrain,
                                   nrow(newdata), trained_bart_BMA$a,trained_bart_BMA$sigma,0,trained_bart_BMA$nu,trained_bart_BMA$lambda,
                                   diff_inital_resids,
                                   all_treated_data,all_control_data,
                                   new_pred_list1)
    }
  }
  
  
  
  
  
  if(is.null(newdata) && length(trained_bart_BMA)==16){
    #ITE_posterior_sum_trees<-gs_chains[[4]]
    ITE_orig_post_sum_trees<-gs_chains[[5]]
    sigma_chains<-gs_chains[[3]]
  }else if(is.null(newdata) && length(trained_bart_BMA)==14){
    #ITE_posterior_sum_trees<-gs_chains[[6]]
    ITE_orig_post_sum_trees<-gs_chains[[7]]
    sigma_chains<-gs_chains[[3]]
    
  }else{
    #ITE_posterior_sum_trees<-gs_chains[[4]]
    ITE_orig_post_sum_trees<-gs_chains[[5]]
    sigma_chains<-gs_chains[[3]]
  } 
  
  
  
  sum_of_tree_BIC<- -0.5*trained_bart_BMA$bic
  weights<-exp(sum_of_tree_BIC-(max(sum_of_tree_BIC)+log(sum(exp(sum_of_tree_BIC-max(sum_of_tree_BIC))))))
  #final_length<-num_iter-burnin
  num_its_to_sample<-round(weights*(num_iter-burnin))
  #final_sigma_chain<-numeric(0)
  
  #final_ITE_chain<-matrix(nrow=0,ncol=ncol(ITE_posterior_sum_trees[[1]]))
  final_ITEorig_chain<-matrix(nrow=0,ncol=ncol(ITE_orig_post_sum_trees[[1]]))
  
  for(i in 1:length(sigma_chains)){
    sample_its<-sample(burnin:num_iter,num_its_to_sample[i])
    #final_sigma_chain<-c(final_sigma_chain,sigma_chains[[i]][sample_its])
    #now do the same for predicted response updates
    #post_ITE_i<-ITE_posterior_sum_trees[[i]]
    post_ITEorig_i<-ITE_orig_post_sum_trees[[i]]
    #final_ITE_chain<-rbind(final_ITE_chain,post_ITE_i[sample_its,])
    final_ITEorig_chain<-rbind(final_ITEorig_chain,post_ITEorig_i[sample_its,])
    
  }
  PI<-apply(final_ITEorig_chain,2,function(x)quantile(x,probs=c(l_quant,0.5,u_quant)))
  meanpreds<-apply(final_ITEorig_chain,2,mean)
  
  CATEs_across_iterations <- apply(final_ITEorig_chain,1,mean)
  CATE_est <- mean(CATEs_across_iterations)
  CATE_PI<-as.matrix(quantile(CATEs_across_iterations,probs=c(l_quant,0.5,u_quant)))
  
  PATE_var<- mean(apply(final_ITEorig_chain,1:2,function(x){(x-CATE_est)^2}))
  PATE_PI <- as.matrix(c(CATE_est + qnorm(l_quant)*sqrt(PATE_var),CATE_est ,  CATE_est + qnorm(u_quant)*sqrt(PATE_var) ))
  
  CATTs_across_iterations <- apply(final_ITEorig_chain[,z_train],1,mean)
  CATT_est <- mean(CATTs_across_iterations)
  CATT_PI<-as.matrix(quantile(CATTs_across_iterations,probs=c(l_quant,0.5,u_quant)))
  
  PATT_var<- mean(apply(final_ITEorig_chain[,z_train],1:2,function(x){(x-CATT_est)^2}))
  PATT_PI <- as.matrix(c(CATT_est + qnorm(l_quant)*sqrt(PATT_var),CATT_est ,  CATT_est + qnorm(u_quant)*sqrt(PATT_var) ))
  
  
  
  CATNTs_across_iterations <- apply(final_ITEorig_chain[,1-z_train],1,mean)
  CATNT_est <- mean(CATNTs_across_iterations)
  CATNT_PI<-as.matrix(quantile(CATNTs_across_iterations,probs=c(l_quant,0.5,u_quant)))
  
  PATNT_var<- mean(apply(final_ITEorig_chain[,1-z_train],1:2,function(x){(x-CATNT_est)^2}))
  PATNT_PI <- as.matrix(c(CATNT_est + qnorm(l_quant)*sqrt(PATNT_var),CATNT_est ,  CATNT_est + qnorm(u_quant)*sqrt(PATNT_var) ))
  
  
  ret<-list()
  ret$ITE_PI<-PI
  ret$ITE_mean <- meanpreds
  ret$CATE_est <- CATE_est
  ret$CATE_PI <- CATE_PI
  ret$CATT_est <- CATT_est
  ret$CATT_PI <- CATT_PI
  ret$CATNT_est <- CATNT_est
  ret$CATNT_PI <- CATNT_PI
  ret$PATE_PI <- PATE_PI
  ret$PATT_PI <- PATT_PI
  ret$PATNT_PI <- PATNT_PI
  ret$fitted.values <- trained_bart_BMA$fitted.values
  
  ret$PIPs <- varIncProb(trained_bart_BMA)
  
  
  class(ret)<-"ITE_intervals.bartBMA"  
  ret
  
}