jointnessBartBMA <- function(object , method=c("DW","LS1","LS2","YQ","YQMa")){
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
    jointness <- log(tl/(bl+br+tl))
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=Inf #this can be an issue for DW when bl==br==0
    return(jointness)
  } else if (method=="LS2"){
    ## Lee-Steel jointness statistic alternative
    jointness <- log(tl/(bl+br))
    jointness <- round(jointness,1)
    diag(jointness) <- NA
    jointness[is.nan(jointness)]=Inf #this can be an issue for DW when bl==br==0
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