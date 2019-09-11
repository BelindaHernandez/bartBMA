#' @title Prediction intervals for bart-bma output obtained using linear algebra to obtain means and variances, and using bisection to find the quantiles of the mixture of t distributions.
#' 
#' @description This function produces prediction intervals for bart-bma output obtained using rmixt function from the package mvnfast.
#' @param min_possible The highest value that the upper quantile of any prediction interval could possibly take. This is required for a bisection root finding algorithm.
#' @param max_possible The lowest value that the upper quantile of any prediction interval could possibly take. This is required for a bisection root finding algorithm.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residuals.
#' @param root_alg_precision The algorithm should obtain approximate bounds that are within the distance root_alg_precision of the true quantile for the chosen average of models.
#' @param spike_tree If equal to 1, then the Spike-and-Tree prior will be used, otherwise the standard BART prior will be used. The number of splitting variables has a beta-binomial prior. The number of terminal nodes has a truncated Poisson prior, and then a uniform prior is placed on the set of valid constructions of trees given the splitting variables and number of terminal nodes.
#' @param lambda_poisson This is a parameter for the Spike-and-Tree prior. It is the parameter for the (truncated and conditional on the number of splitting variables) Poisson prior on the number of terminal nodes.
#' @export 
#' @return The output is a list of length one. The one element in this list is a vector of prediction intervals???

bartBMA_with_ITEs_exact_par <-function(l_quant,u_quant,newdata=NULL,update_resids=1,
                                       num_cores=1,root_alg_precision=0.00001,
                                       x_covariates,z_train ,y_train,
                                       a1=3,nu1=3,sigquant1=0.9,c1=1000,
                                       pen1=12,num_cp1=20,x.test1=matrix(0.0,0,0),
                                       num_rounds1=5,alpha1=0.95,beta1=2,split_rule_node1=0,
                                       gridpoint1=0,maxOWsize1=100,num_splits1=5,
                                       gridsize1=10,zero_split1=1,only_max_num_trees1=1,
                                       min_num_obs_for_split1=2, min_num_obs_after_split1=2,
                                       exact_residuals1=1,spike_tree1=0,lambda_poisson1=10,less_greedy1=0){
  
  
  x_train <- cbind(z_train,x_covariates)
  
  object <- bartBMA(x.train= x_train,y.train =y_train ,
                    a=a1,nu=nu1,sigquant=sigquant1,c=c1,
                    pen=pen1,num_cp=num_cp1,x.test=x.test1,
                    num_rounds=num_rounds1,alpha=alpha1,beta=beta1,split_rule_node=split_rule_node1,
                    gridpoint=gridpoint1,maxOWsize=maxOWsize1,num_splits=num_splits1,
                    gridsize=gridsize1,zero_split=zero_split1,only_max_num_trees=only_max_num_trees1,
                    min_num_obs_for_split=min_num_obs_for_split1, min_num_obs_after_split=min_num_obs_after_split1,
                    exact_residuals=exact_residuals1,spike_tree=spike_tree1,lambda_poisson=lambda_poisson1,less_greedy=less_greedy1)
  
  
  
  if(is.null(newdata) && length(object)==16){
    #if test data specified separately
    ret<-pred_ints_ITE_outsamp_par(object$sumoftrees,object$obs_to_termNodesMatrix,
                                   object$response,object$bic,#min_possible, max_possible,
                                   object$nrowTrain,
                                   nrow(object$test_data),object$a,object$sigma,0,object$nu,
                                   object$lambda,#diff_inital_resids,
                                   object$test_data,l_quant,u_quant,num_cores,
                                   root_alg_precision,x_covariates
    )
    
    
    
  }else{if(is.null(newdata) && length(object)==14){
    #else return Pred Ints for training data
    ret <- pred_ints_ITE_insamp_par(object$sumoftrees,
                                    object$obs_to_termNodesMatrix,
                                    object$response,object$bic,#min_possible, max_possible,
                                    object$nrowTrain,#nrow(object$test_data),
                                    object$a,object$sigma,0,object$nu,
                                    object$lambda,#diff_inital_resids,object$test_data,
                                    l_quant,u_quant,
                                    num_cores,
                                    root_alg_precision,x_covariates
    )
    
  }else{
    #if test data included in call to object
    ret<-pred_ints_ITE_outsamp_par(object$sumoftrees,object$obs_to_termNodesMatrix,object$response,object$bic,#min_possible, max_possible,
                                   object$nrowTrain,
                                   nrow(newdata), object$a,object$sigma,0,object$nu,
                                   object$lambda,#diff_inital_resids,
                                   newdata,l_quant,u_quant,num_cores,
                                   root_alg_precision,x_covariates
    )
    
  }}
  
  #PI<-apply(draws_from_mixture,2,function(x)quantile(x,probs=c(l_quant,0.5,u_quant)))
  
  
  
  #each row is a vector drawn from the mixture distribution
  
  
  names(ret)<-c("ITE_intervals",
                "ITE_estimates",
                "CATE_Interval",
                "CATE_estimate")
  
  
  ret
}