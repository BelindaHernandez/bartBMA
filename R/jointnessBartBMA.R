#' Jointness measure for all pairs of variables in bartBMA output
#'
#' `jointnessBartBMA()` returns a matrix giving a jointness measure for all variable pairs. The variables are those in the covariate matrix used to create the bartBMA object, and are in the same order as the columns of the original covariate matrix.
#'
#' This function obtains jointness measures for a bartBMA object. Jointness
#' is based on the joint posterior distribution of variables over the model
#' space. Positive jointness implies that covariates are complements,
#' representing distinct but mutually reinforcing effects. Negative jointness
#' implies that explanatory variables are substitutes, and capture similar
#' underlying effects \insertCite{doppelhofer2009jointnessa}{bartBMA}. \cr
#' \insertCite{doppelhofer2009jointnessa;textual}{bartBMA} introduce a jointness
#' measure defined as follows 
#' \deqn{\ln \left[ \frac{p(i \cap l | y)p(\bar{i} \cap \bar{l}|y)}{p(i\cap \bar{l}|y)p(\bar{i}\cap l|y)}\right]}
#' where \eqn{p(i \cap l | y)} denotes the posterior joint inclusion probability of variables \eqn{i} and \eqn{l}
#' and \eqn{p(\bar{i} \cap \bar{l}|y)}denotes the posterior joint exclusion probability of variables \eqn{i} and \eqn{l},
#' and so on.\cr
#' \insertCite{ley2007jointness;textual}{bartBMA} introduce the following jointness measures:\cr
#' This function refers to the following measure as LS1:
#' \deqn{ \left[ \frac{p(i \cap l | y)}{p(i \cap \bar{l}|y)+p(\bar{i}\cap l|y)+p(i \cap l | y)} }
#' and to the following measure as LS2:
#' \deqn{ \left[ \frac{p(i \cap l | y)}{p(i \cap \bar{l}|y)+p(\bar{i}\cap l|y)}}
#' LS1 and LS2 are always non-negative, and are defined as long as at least
#' one variable has nonzero inclusion probability.\cr
#' \insertCite{hofmarcher2018bivariate;textual}{bartBMA} and 
#' \insertCite{cuaresma2015comprehensive;textual}{bartBMA} propose the use of a transformation
#' of DW to the interval \eqn{[-1,1]}, which is Yule's Q
#' \deqn{\frac{p( i \cap l|y)p(\bar{i} \cap \bar{l}|y)-p(i \cap \bar{l}|y)p(\bar{i} \cap l|y)}{p( i \cap l|y)p(\bar{i} \cap \bar{l}|y)+p(i \cap \bar{l}|y)p(\bar{i} \cap l|y)}}
#' \insertCite{hofmarcher2018bivariate;textual}{bartBMA} note that YQ and DW are different
#' to LS1 and LS2 in that YQ and DW are non null invariant, i.e. they take
#' account of when both variables are excluded. This corresponds to the observation by
#' \insertCite{doppelhofer2009jointnessb;textual}{bartBMA} that jointness can manifest itself in both the inclusion and exclusion margin of the joint posterior distribution.\cr
#' \insertCite{hofmarcher2018bivariate;textual}{bartBMA} argue in favour of non null-variant measures,
#' which allow pairs of variables to be categorised into substitututes and complements.\cr
#' \insertCite{hofmarcher2018bivariate;textual}{bartBMA} also propose a new measure, Modified Yule's Q (YQM),
#' which is always defined, and satisfies some additional desirable proporties.
#' However, Modified Yule's Q is only applicable in the context of BMA implemented by
#' Markov Chain Monte Carlo Model Composition, and we include here a measure that
#' makes a similar adjustment, YQMa, but which can be implemented for BMA more generally\;
#' First, let \eqn{\tilde{p}(i \cap l|y) = p(i \cap l |y) +\frac{\frac{1}{2}}{N}}, 
#' where \eqn{N} is the number of iterations for YQM, but is some arbitraty large number
#' for YQMa. We set this number to 5,000,000. Then YQMa can be defined as
#' \deqn{\frac{N^2 \left[ \tilde{p}( i \cap l|y)\tilde{p}(\bar{i} \cap \bar{l}|y)-\tilde{p}(i \cap \bar{l}|y)\tilde{p}(\bar{i} \cap l|y) \right]}{N^2 \left[ \tilde{p}( i \cap l|y)\tilde{p}(\bar{i} \cap \bar{l}|y)+\tilde{p}(i \cap \bar{l}|y)\tilde{p}(\bar{i} \cap l|y) \right]-\frac{1}{2}} }
#' which can also be written as
#' \deqn{\frac{\left[ \tilde{p}( i \cap l|y)\tilde{p}(\bar{i} \cap \bar{l}|y)-\tilde{p}(i \cap \bar{l}|y)\tilde{p}(\bar{i} \cap l|y) \right]}{ \left[ \tilde{p}( i \cap l|y)\tilde{p}(\bar{i} \cap \bar{l}|y)+\tilde{p}(i \cap \bar{l}|y)\tilde{p}(\bar{i} \cap l|y) \right]-\frac{1}{2N^2}}}
#' YQMa is always defined, whereas YQM and DW are undefined when one of the variables is always
#' included or one of the variables is always excluded. \insertCite{doppelhofer2009jointnessa;textual}{bartBMA}
#' introduce a modified version of DW, denoted here as DWM, such that if only
#' one of the variables is always included or excluded, then DWM is the log of the
#' posterior odds of including the other variable. For example, if \eqn{l} is always included
#' or always excluded, then the jointness between \eqn{i} and \eqn{l} is \eqn{\ln \left(  \frac{p(i | y)}{1-p(i|y)} \right)} . \cr
#' \insertCite{man2018critical;textual}{bartBMA} proposes an alternative measure, which simply takes
#' the average of one of LS1 or LS2, and DW. The option provided by this package is the average of
#' LS1 and DW, and is equal to LS1 if DW is undefined.\cr
#' Another null invariant measure, proposed by \insertCite{strachan2009comment;textual}{bartBMA}
#' is:
#' \deqn{p(i|y)p(l|y)\ln \left(\frac{p(i \cap l |y)}{p(i \cap \bar{l} |y)p( \bar{i} \cap l |y)} \right) }
#' 
#'
#' @param object A bartBMA object
#' @param method As detailed in the description, there are a number of definitions of jointness
#' \itemize{
#' \item LS1 and LS2 are methods defined by \insertCite{ley2007jointness;textual}{bartBMA}.
#' \item DW is the method defined by \insertCite{doppelhofer2009jointnessa;textual}{bartBMA}.
#' \item DWM is a modification of DW suggested in footnote 12 of the original paper by \insertCite{doppelhofer2009jointnessa;textual}{bartBMA}.
#' \item Str is defined by \insertCite{strachan2009comment;textual}{bartBMA}
#' \item YQ is modified Yule's Q, introduced by \insertCite{hofmarcher2018bivariate;textual}{bartBMA} and \insertCite{cuaresma2015comprehensive;textual}{bartBMA}.
#' \item YQMa is similar to modified Yule's Q introduced by \insertCite{hofmarcher2018bivariate;textual}{bartBMA}. YQMa is always defined.
#' \item Man is the average of LS2 and DW, suggested by \insertCite{man2018critical;textual}{bartBMA}
#' }
#' @return The output is a symmetric matrix, with rows and columns in the same order as the columns of the
#' covariate matrix used to create the bartBMA object. Each off-digonal element is the jointness between a pair of variables.
#'
#' @examples
#' \dontrun{
#' N <- 100
#' set.seed(100)
#' x1 <- runif(N)
#' x2 <- runif(N)
#' x3 <- runif(N)
#' x4 <- runif(N)
#' x5 <- runif(N)
#' x6 <- runif(N)
#' x7 <- runif(N)
#' x8 <- runif(N)
#' x9 <- runif(N)
#' x10 <- runif(N)
#' epsilon <- rnorm(N)
#' xcov <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
#' y <- sin(pi*x1*x2) + 20*(x3-0.5)^2+10*x4+5*x5+epsilon
#' bart_bma_example <- bartBMA(x.train = xcov, y.train = y)
#' jointnessBartBMA(bart_bma_example, method ="YQMa")
#' }
#' @author Eoghan O'Neill
#' @references
#' \insertAllCited{}
#' @export
#' @importFrom Rdpack reprompt 
#' @importFrom Rcpp evalCpp

jointnessBartBMA <- function(object , method=c("DW","LS1","LS2","YQ","YQMa","DWM","Man","Str")){
  if (class(object) != "bartBMA") { #user check
    if (is.data.frame(object) | is.matrix(object)) {
      stop("need  a bartBMA object. Run bartBMA on the data first")
    } else stop("need  a bartBMA object.")
  }
  method=method[[1]] # User check
  
  imp_vars2=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees)
  pmps <- imp_vars2[[1]] #Posterior Model Probabilities (BIC weights)
  tl <- crossprod(((imp_vars2[[3]]>0)*pmps),((imp_vars2[[3]]>0))) #PMP where both x&y are included
  tr <- crossprod(((imp_vars2[[3]]==0)*pmps),((imp_vars2[[3]]==0))) #PMP where both x&y are included
  bl <- crossprod(((imp_vars2[[3]]==0)*pmps),((imp_vars2[[3]]>0))) #PMP where both x&y are included
  br <- crossprod(((imp_vars2[[3]]>0)*pmps),((imp_vars2[[3]]==0))) #PMP where both x is included and y is excluded
  
  if (method=="LS1"){
    ## Lee-Steel jointness statistic
    jointness <- (tl/(bl+br+tl))
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #undefined if both variables always excluded
    return(jointness)
  } else if (method=="LS2"){
    ## Lee-Steel jointness statistic alternative
    jointness <- (tl/(bl+br))
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #undefined if both variables always excluded
    return(jointness)
  } else if (method=="Man"){
    ## Lee-Steel jointness statistic alternative
    jointness1 <- (tl/(bl+br+tl))
    diag(jointness1) <- NA
    jointness1[is.nan(jointness1)]=NA #undefined if both variables always excluded
    jointness2 <- log((tl*tr)/(bl*br))
    diag(jointness2) <- NA
    jointness2[is.nan(jointness2)]=NA #this can be an issue for DW when bl==0 or br==0
    jointness <- apply(array(c(jointness1, jointness2),dim=c(nrow(jointness1),ncol(jointness1),2)),1:2,mean, na.rm=TRUE)
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #this can be an issue for DW when bl==0 or br==0
    return(jointness)
  } else if (method=="Str"){
    incProbs <- varIncProb(object) #obtain single variable inclusion probabilities
    jointness <- t(t(incProbs*((tl)/(bl*br)))*incProbs)
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #this can be an issue for YQ when tl*tr==0 and bl*br==0 # note, numerator can be negative, so can't make plus inf as for DW
    return(jointness)
  } else if (method=="DWM"){
    incProbs <- varIncProb(object) #obtain single variable inclusion probabilities
    incProbMatrixi <- replicate(length(incProbs),incProbs)
    incProbMatrixj <- t(incProbMatrixi)
    jointness <- ifelse((incProbMatrixi== 0|incProbMatrixi==1)&(incProbMatrixj== 0|incProbMatrixj==1),NA,
                        ifelse((incProbMatrixi== 0|incProbMatrixi==1),incProbMatrixj/(1-incProbMatrixj),
                               ifelse((incProbMatrixj== 0|incProbMatrixj==1) , incProbMatrixi/(1-incProbMatrixi),
                                      log((tl*tr)/(bl*br))) ) )
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #this can be an issue for YQ when tl*tr==0 and bl*br==0 # note, numerator can be negative, so can't make plus inf as for DW
    return(jointness)
  } else if (method=="YQ"){
    ## Yule's Q jointness statistic alternative
    jointness <- (tl*tr-bl*br)/(tl*tr+bl*br)
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #this can be an issue for YQ when tl*tr==0 and bl*br==0 # note, numerator can be negative, so can't make plus inf as for DW
    return(jointness)
  } else if (method=="YQMa"){
    #arbitrary number to make YQM feasible
    #For MC^3 methods this would be the number of iterations
    fake_iternum <- 5000000
    
    #adjustd joint inclusion probs for YQM
    tla <- tl+(0.5/fake_iternum)
    tra <- tr+(0.5/fake_iternum)
    bla <- bl+(0.5/fake_iternum)
    bra <- br+(0.5/fake_iternum)
    ## Yule's Q jointness statistic
    jointness2 <- (tla*tra-bla*bra)/(tla*tra+bla*bra - (1/(2*(fake_iternum^2))))
    jointness2 <- round(jointness2,1)
    diag(jointness2) <- "."
    jointness2 <- as.data.frame(jointness2)
    return(jointness2)
  } else {
    ## Doppelhofer-Weeks jointness statistic
    jointness <- log((tl*tr)/(bl*br))
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=NA #this can be an issue for DW when bl==0 or br==0
    return(jointness)
  }
}