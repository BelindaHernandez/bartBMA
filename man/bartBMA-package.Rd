\name{bartBMA-package}
\alias{bartBMA}
\docType{package}
\title{
Bayesian Additive Regression Trees using Bayesian Model Averaging
}
\description{
\packageDescription{bartBMA}
}
\usage{
bartBMA(x.train,y.train,a=3,nu=3,sigquant=0.9,c=1000,
                     pen=12,num_cp=20,x.test=matrix(0.0,0,0),num_rounds=5, alpha=0.95, beta=1,split_rule_node=0,gridpoint=0)
}
\arguments{
   \item{x.train}{
   Numeric matrix of explanatory variables for training (in sample) data,
   with rows corresponding to observations and columns to variables.\cr
   }
   \item{y.train}{
   Response variable for training (in sample) observations. Response can either be numeric or binary.\cr
   If y is numeric then a numeric response model will be fit.\cr
   If y is binary (0,1) then a probit link will be used. \cr
   }
   \item{a}{
   Hyperparameter on terminal node mean variance \eqn{\mu_{ij}}. Here the same terminal node prior as Chipman, George and McCulloch (1998) is used where by \eqn{\mu_{ij} \sim N(0,\sigma^2/a)}. Default value is 3.  
   }
   \item{nu}{
   Value of hyperparameter for residual precision \eqn{\tau} which has a \eqn{\Gamma(\nu/2, \nu\lambda/2)} prior.Default is 3.
   }
   \item{sigquant}{
   The quantile of the prior that the rough estimate of the prior standard deviation is placed at.
   Used to set the prior on hyperparameter \eqn{\lambda} where the residual precision \eqn{\tau} has a \eqn{\Gamma(\nu/2, \nu\lambda/2)} prior. Default value is 0.9.
   }
   \item{c}{
   Determines the size of Occam's window. Any sum of tree model which fits the following critera are saved to the final model otherwise they are deleted: \eqn{BIC_\ell-argmin_\ell(BIC_\ell)\leq log(c)}. Default value is 1000.
   }
   \item{pen}{
   If \code{gridpoint}=0, \code{pen} determines the penalty parameter for adding additional change points in the \code{PELT} algorithm. If \eqn{n <200} a good default value is \eqn{10 \log n}. If \eqn{n <200} it is recommended to use the gridpoint search to find the best set of split rules. \cr  
   If \code{gridpoint}=1, \code{pen} determines the number of grid points to be searched within the range of each variable in dataset \code{x.train}.\cr Default value 12.
   }
   \item{num_cp}{
   The \eqn{\%} of best split rules which are saved to memory and used to build trees. Default value is 20.
   }
   \item{x.test}{
   Numeric matrix of explanatory variables for training (out of sample) data,
   with rows corresponding to observations and columns to variables.\cr 
   Should have the same structure and variables as x.train.
   }
   \item{num_rounds}{
   Maximum length of each sum of tree model. Default is 5 meaning each sum of tree model can contain a maximum of 5 decision trees.
   }
   \item{alpha}{
   Hyperparameter on prior probability that a given node in a tree is non-terminal. Default value \eqn{\alpha}=0.95.
   }
   \item{beta}{
   Hyperparameter on prior probability that a given node in a tree is non-terminal. Default value \eqn{\alpha}=1.
   }
   \item{split_rule_node}{
   Binary indicator to determine whether the set of best split rules should be searched for before growing trees or should be updated at each internal node
   of each tree.\cr
   If split_rule_node=0 the best set of split rules will be searched before the each tree is grown.\cr
   If split_rule_node=1 the best set of split rules will be updated at each internal node as each tree is being grown.\cr
   }
   \item{gridpoint}{
   Binary indicator to determine method of searching for best split rules when growing trees.\cr
   If gridpoint=0 then PELT (Pruned Exact Linear time) changepoint detection algorithm will be used.\cr
   If gridpoint=1 then a grid search will be performed using pen equally spaced grid points within the range of the variable being considered. \cr
   }   
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bartBMA}
\packageIndices{bartBMA}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{bartBMA}

Maintainer: \packageMaintainer{bartBMA}
}
\references{
Hernandez B, Raftery AE, Pennington SR, Parnell AC,BART-BMA Bayesian Additive Regression Trees using Bayesian Model Averaging, Journal of Statistics and Computing
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }

