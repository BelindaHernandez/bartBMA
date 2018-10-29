#create main function call and bartBMA class

bartBMA<-function(x,...)UseMethod("bartBMA")

bartBMA.default<-function(x.train,y.train,a=3,nu=3,sigquant=0.9,c=1000,
                          pen=12,num_cp=20,x.test=matrix(0.0,0,0),num_rounds=5,alpha=0.95,beta=1,split_rule_node=0,gridpoint=0,maxOWsize=100){
  binary=FALSE
  start_mean=0
  start_sd=1
  mu=0
  sigma_mu=0; 
  sigma=sd(y.train)/(max(y.train)-min(y.train))
  qchi = qchisq(1.0-sigquant,nu,1,0);
  lambda = (sigma*sigma*qchi)/nu;
  if(is.factor(y.train)) {
    # if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
    binary = TRUE
    #  y.train = as.numeric(y.train)-1
    stop("Response must be a numeric vector")
  } else {
    if((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0)) {
      cat('NOTE: assumming numeric response is binary\n')
      binary = TRUE
      #stop("Response must be a numeric vector")
    }
  }
  
  if(is.vector(x.train) | is.factor(x.train)| is.data.frame(x.train)) x.train = as.matrix(x.train)
  if(is.vector(x.test) | is.factor(x.test)| is.data.frame(x.train)) x.test = as.matrix(x.test)
  
  if(is.matrix(x.train)) {
    if(nrow(x.test)) {
      if(!is.matrix(x.test)) stop('x.test must be a matrix')
    } 
  }
  #check input arguments:
  # if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
  #if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
  if((!is.matrix(x.train))) stop("argument x.train must be a double matrix")
  if((!is.matrix(x.test)) ) stop("argument x.test must be a double matrix")
  if(!binary) {
    if((!is.vector(y.train))) stop("argument y.train must be a double vector")
  }
  if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
  if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
  #if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
  #if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
  #if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
  if(c<1)stop("Value of Occam's Window has to be greater than 0."); 
  if(num_cp<0 || num_cp>100)stop("Value of num_cp should be a value between 1 and 100."); 
  
  bartBMA_call=BART_BMA_sumLikelihood(x.train,y.train,start_mean,start_sd,a,mu,nu,lambda,c,sigma_mu,
                                      pen,num_cp,x.test,num_rounds,alpha,beta,split_rule_node,gridpoint,maxOWsize)
  
  if(length(bartBMA_call)==6){
    #length of bartBMA_call is 6 if test data was included in the call
    names(bartBMA_call)<-c("fitted.values","sumoftrees","obs_to_termNodesMatrix","bic","test.preds","sum_residuals")
    bartBMA_call[[6]]<-bartBMA_call[[6]][[length(bartBMA_call[[6]])]]
    bartBMA_call$test_data<-x.test
  }else{
    names(bartBMA_call)<-c("fitted.values","sumoftrees","obs_to_termNodesMatrix","bic","sum_residuals")
    bartBMA_call[[5]]<-bartBMA_call[[5]][[length(bartBMA_call[[5]])]]
  }
  
  bartBMA_call$numvars<-ncol(x.train)
  bartBMA_call$call<-match.call()
  bartBMA_call[[2]]<-bartBMA_call[[2]][[length(bartBMA_call[[2]])]]
  bartBMA_call[[3]]<-bartBMA_call[[3]][[length(bartBMA_call[[3]])]]
  bartBMA_call$y_minmax<-range(y.train)
  bartBMA_call$response<-y.train
  bartBMA_call$nrowTrain<-nrow(x.train)
  bartBMA_call$sigma<-sigma
  bartBMA_call$a<-a
  bartBMA_call$nu<-nu
  bartBMA_call$lambda<-lambda
  
  class(bartBMA_call)<-"bartBMA"
  bartBMA_call
}