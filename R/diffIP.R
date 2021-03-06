#' @title diffIP
#' @param x The RNADMethyl data list
#' @param Covariates Matrix of covariates to include in the PoissonGamma test. Default is NULL.
#' @param p.adjBy The method used to adjust for the p value. Default is 'fdr'. Can be "bonferroni","BY","fdr".
#' @param exclude Names of samples to exclude (as outlier) in the test. Default is NULL
#' @param maxPsi The cutoff value for random effect parameter Psi estimation
#' @param mode Run test on "bin" mode or "peak" mode.
#' @return x The RNADMethyl data list with estimators
#' @export
diffIP <- function(
  x,
  Covariates = NULL, # covariates
  plotPvalue = TRUE,
  p.adjBy = "fdr", # the method to calculate adjusted p_value
  exclude = NULL,
  mode,
  maxPsi = 100
){
  #### start run PoissonGamma test
  if( is.null(Covariates) ){ # run simple PoissonGamma

    if(is.null(exclude) ){ ## running default
      if(mode == "bin"){
        allY <- x$ip_adjExpr_filtered
      }else if(mode == "peak"){
        allY <- x$jointPeak_adjExpr
      }else{
        stop("Must choose the mode, either bin mode or peak mode.")
      }

      X <- x$X
    }else if( all( exclude %in% x$samplenames) ){
      exc.id <- match(exclude,x$samplenames)
      if(mode == "bin"){
        allY <- x$ip_adjExpr_filtered[,-c(exc.id)]
      }else if(mode == "peak"){
        allY <- x$jointPeak_adjExpr[,-c(exc.id)]
      }else{
        stop("Must choose the mode, either bin mode or peak mode.")
      }
      X <- x$X[-c(exc.id)]
      cat(paste0("Sample ",exclude, " has been removed from test...\n"))
    }else{
      stop("Samples to exclude is not in samplenames...")
    }

    tmp <- rafalib::as.fumeric(X) - 1
    names(tmp) <- X
    cat("The categorical variable has been converted:\n")
    print(tmp)
    X <- rafalib::as.fumeric(X) - 1 # convert categorical variable into numerical variable.


    psi <- 10

    cat("running PoissonGamma test at single beta mode...\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all.est <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ])
      model1 <- glm(Y ~ X, family = poisson(link = 'log'))
      coef <- model1$coefficients
      mu2 <- coef[1]
      beta <- coef[2]
      est <- try(unlist(PoissionGamma(Y, X, beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all.est <- rbind(all.est, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }

    rownames(all.est) <- rownames(allY)[all.id]
    if(plotPvalue){ # plot the distribution of p values
      hist(all.est[,"p_value"],main = "P Values", xlab = "p_value")
    }

    ##calcualte adjusted p value
    padj <- p.adjust(all.est[,"p_value"],method = p.adjBy )
    all.est <- cbind( all.est, padj )
    if(mode == "bin"){
      x <- c(x,list('all.est'=all.est))
    }else if(mode == "peak"){
      x <- c(x,list('all.est.peak'=all.est))
    }

    return(x)

  } else if( !is.null(Covariates) & dim( cbind(X,Covariates) )[1] == length(X)  ){ # include covariates

    if(is.null(exclude) ){ ## running default
      if(mode == "bin"){
        allY <- x$ip_adjExpr_filtered
      }else if(mode == "peak"){
        allY <- x$jointPeak_adjExpr
      }else{
        stop("Must choose the mode, either bin mode or peak mode.")
      }

      X <- x$X
    }else if( all( exclude %in% x$samplenames) ){
      exc.id <- match(exclude,x$samplenames)
      if(mode == "bin"){
        allY <- x$ip_adjExpr_filtered[,-c(exc.id)]
      }else if(mode == "peak"){
        allY <- x$jointPeak_adjExpr[,-c(exc.id)]
      }else{
        stop("Must choose the mode, either bin mode or peak mode.")
      }
      X <- x$X[-c(exc.id)]
      Covariates <- Covariates[-c(exc.id),]
      cat(paste0("Sample ",exclude, " has been removed from test...\n"))
    }else{
      stop("Samples to exclude is not in samplenames...")
    }

    tmp <- rafalib::as.fumeric(X) - 1
    names(tmp) <- X
    cat("The categorical variable has been converted:\n")
    print(tmp)
    X1 <- rafalib::as.fumeric(X) - 1 # convert categorical variable into numerical variable.

    ## organize covariates
    if(!is.numeric(Covariates)){stop("Please convert covariates into numerical variables")}

    X.all <- cbind(X1,Covariates) # new design matrix
    colnames(X.all) <- paste("X",1:ncol(X.all),sep = "")

    #design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[", rownames(X.all), sep = "",collapse = " + ")) )
    design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[,", 1:ncol(X.all),"]", sep = "",collapse = " + ")) )
    ## Run multi-beta PoissonGamma
    psi <- 10

    cat("running PoissonGamma test at multi-beta mode...\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all1 <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ] )
      aa <- unlist(summary( lm( design.multiBeta ) )$coefficients[, 1])
      mu2 <- aa[1]
      beta <- aa[2:(ncol(X.all)+1 )]
      est <- try(unlist(PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = 0.25, steps = 10, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all1 <- rbind(all1, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }

    rownames(all1) <- rownames(allY)[all.id] ## assign window names to test statistics
    if(plotPvalue){ # plot the distribution of p values
      hist(all1[,"p_value3"],main = "P Values", xlab = "p_value")
    }
    ##calcualte adjusted p value
    padj <- p.adjust(all1[,"p_value3"],method = p.adjBy )
    all1 <- cbind( all1, padj )
    if(mode == "bin"){
      x <- c(x,list('all.est'=all1))
    }else if(mode == "peak"){
      x <- c(x,list('all.est.peak'=all1))
    }

    return(x)

  } else{ ## if Covariate is not null and not the same dimension of X, throw a error
    stop("Covariates need to be the same dimension as X")
  }

}
