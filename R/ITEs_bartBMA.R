#' @title ITE Predictions (in-sample) using bartBMA and the method described by Hill (2011)
#' 
#' @description This function produces ITE Predictions (in-sample) using bartBMA and the method described by Hill (2011).
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
#' @export 
#' @return A vector of Individual Treatment Effect Estimates.

ITEs_bartBMA<-function(x_covariates,z_train ,y_train,
                       a1=3,nu1=3,sigquant1=0.9,c1=1000,
                       pen1=12,num_cp1=20,x.test1=matrix(0.0,0,0),
                       num_rounds1=5,alpha1=0.95,beta1=2,split_rule_node1=0,
                       gridpoint1=0,maxOWsize1=100,num_splits1=5,gridsize1=10,zero_split1=1,only_max_num_trees1=1,
                       min_num_obs_for_split1=2, min_num_obs_after_split1=2){
  
  x_train <- cbind(z_train,x_covariates)
  
  trained_bart_BMA <- bartBMA(x.train= x_train,y.train =y_train ,
    a=a1,nu=nu1,sigquant=sigquant1,c=c1,
    pen=pen1,num_cp=num_cp1,x.test=x.test1,
    num_rounds=num_rounds1,alpha=alpha1,beta=beta1,split_rule_node=split_rule_node1,
    gridpoint=gridpoint1,maxOWsize=maxOWsize1,num_splits=num_splits1,gridsize=gridsize1,zero_split=zero_split1,only_max_num_trees=only_max_num_trees1,
    min_num_obs_for_split=min_num_obs_for_split1, min_num_obs_after_split=min_num_obs_after_split1)
  
  if(nrow(x.test1)==0){
    all_treated_data <- cbind(rep(1,nrow(x_covariates)), x_covariates)
    all_control_data <- cbind(rep(0,nrow(x_covariates)), x_covariates)
    
    preds_treated<-get_BART_BMA_test_predictions(all_treated_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
    preds_control<-get_BART_BMA_test_predictions(all_control_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
  }else{
    all_treated_data <- cbind(rep(1,nrow(x.test1)), x.test1)
    all_control_data <- cbind(rep(0,nrow(x.test1)), x.test1)
    
    preds_treated<-get_BART_BMA_test_predictions(all_treated_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
    preds_control<-get_BART_BMA_test_predictions(all_control_data,trained_bart_BMA$bic,trained_bart_BMA$sumoftrees,trained_bart_BMA$y_minmax)
    
  }
  
  ITE_ests<-preds_treated[[1]]-preds_control[[1]]
  class(ITE_ests)<-"ITE_ests.bartBMA"
  ITE_ests
}