

#' Extended Bayesian Information Criterion for \code{glmnet} objects
#'
#' @param x The \code{glmnet} model.
#' @param gamma The gamma parameter. This should be a value between 0 and 1.
#' At a gamma of 0, EBIC is equivalent to BIC, while at higher values, it is
#' more strict on the number of predictors modelled.
#'
#' @returns The EBIC of the model.
#' @export
#'
#' @examples
#' 
#' 
EBIC.glmnet <- function(x, gamma){
  tLL <- -deviance(x)
  k <- x$df
  n <- x$nobs
  p <- nrow(x$beta)
  EBIC <- -tLL + log(n)*k + 2*gamma*k*log(p)
  return(EBIC)
}

###################################
# Helper functions for regressors #
###################################

getStateMask <- function(currstate, pat, sampletype, direction){
  pat$idx <- pat$idx * direction
  if(sampletype == "all"){
    include_nonvisitors <- T
  }else if(sampletype == "visitors"){
    include_nonvisitors <- F
  }
  mask <- pat |> 
    group_by(.data$sequence) |> 
    mutate(limited=any(.data$state == currstate)) |>
    transmute(mask= ifelse(.data$limited, 
                           .data$idx <= .data$idx[.data$state == currstate], 
                           include_nonvisitors))
  return(mask$mask)
}

# Assemble betas and ebics into matrices
assembleCoefficients <- function(statefits, extrapreds=NULL){
  # Extract ebics
  ebicmat <- sapply(statefits, \(x)x$ebics)
  
  # Add zeroes for each fit's own state and for start
  statenames <- names(statefits)
  betalist <- statefits |> lapply(\(x)x[["betas"]])
  for(statename in statenames){
    # Add empty row for own state
    rowid <- nrow(betalist[[statename]])+1
    betalist[[statename]] %<>% rbind(0)
    rownames(betalist[[statename]])[rowid] <- statename
    
    # Reorder
    betalist[[statename]] <- 
      betalist[[statename]][c("(Intercept)", extrapreds, statenames),]
  }
  
  # Extract betas by gamma
  gammanames <- rownames(ebicmat)
  gamma_coeflist <-
    vector(length=length(gammanames), mode="list") |> 
    setNames(gammanames)
  for(currgammaname in gammanames){
    currbetas <- betalist |> sapply(\(x)x[, currgammaname])
    colnames(currbetas) <- names(betalist)
    gamma_coeflist[[currgammaname]] <- currbetas
  }
  
  # Return output
  out <- list(gammacoefs=gamma_coeflist, ebicmat=ebicmat)
  return(out)
}

extract_coefs <- function(fit, gammas){
  iterminebics <- iterebics <- numeric(length(gammas))
  for(gammidx in seq_along(gammas)){
    ebicvec <- EBIC.glmnet(fit, gamma=gammas[gammidx])
    iterminebics[gammidx] <- which.min(ebicvec)
    iterebics[gammidx] <- min(ebicvec)
  }
  betamat <- coef(fit)[,iterminebics]
  colnames(betamat) <- names(iterebics) <- paste0("gamma", gammas)
  
  out <- list(betas=betamat, ebics=iterebics)
  return(out)
}

################################
# TRANSITION NETWORK FUNCTIONS #
################################

#' Model a transition network from sequence data
#'
#' @param data A data.frame with 3 columns: sequence (character), 
#' state (character), and date (numeric).
#' @param predictors Which states should be included in the network?
#' @param decaytype Once a state has occurred, how should its predictor look?
#' 
#' * If "flat", the predictor is 1 after the state has occurred until it 
#' no longer satisfies the condition set by the parameter argument.
#' 
#' * If "step", the predictor is 1 immediately after the state has occurred,
#' and decays to 0 until it no longer satisfies the condition set by 
#' the parameter argument.
#' 
#' * If "accrual", the predictor is 1 after the state has occurred 
#' until sequence end.
#' @param parameter This determines the largest temporal distance from 
#' the occurrence of a state that is transformed into a nonzero value 
#' as a predictor. The exact behavior of this argument is determined by
#' the \code{decaytype} argument. See Details for clarification.
#' @param predtype What type of predictor should be used? See Details.
#' @param direction In which direction should predictions be made?
#' 1 means current states predict future states, -1 means current states predict
#' past states.
#' @param sampletype When the occurrence of a state is predicted, should every 
#' sequence be included ("all") or only those that contain the state ("visitors")?
#' @param gammas Which EBIC gamma values should be used for model selection?
#' Higher values mean less edges are retained in the network.
#' @param alpha This is the elasticnet mixing parameter.
#' 1 (default) means lasso regularization, 0 means ridge regularization, and 
#' values in-between mix the two. See [glmnet::glmnet()].
#' @param force.positive Should model beta values be restricted to 
#' the positive range?
#' @param ncores Number of models to run in parallel. 
#' Be careful not to overload system memory.
#' @param verbose Produce verbose output in console?
#' 
#' @details
#' # Rationale
#' 
#' This is a method for analyzing sequences of states. The first occurrence 
#' of each state within a sequence (its _discovery_, for short) 
#' is predicted with the recently preceding first or most recent occurrences 
#' of other states, using lasso-regularized logistic regressions. 
#' The \code{preds} argument sets which states are predicted and 
#' used as predictors. All other states are kept in the data, 
#' but do not contribute to prediction nor are predicted.
#' 
#' # Dependent variable
#' 
#' The level of analysis (and the actual data analyzed) is 
#' the discoveries of all states. That is, the imputed sequences 
#' in \code{data} are transformed such that each row represents a discovery, 
#' not a mere occurrence of a state.
#' For each predicted state, the dependent variable is coded as 0 
#' if the discovery on the current row is not the predicted state, 
#' and 1 if it is that state; prediction within a sequence continues until 
#' this state or the end of the sequence is reached. 
#' Hence, if a state is discovered before the end of the sequence, 
#' all discoveries after it are ignored in the prediction of that state, 
#' since there can only be one discovery of a state in a sequence,
#' and therefore it would not make sense to keep predicting a discovery
#' after it has already occurred.
#' 
#' # Predictors
#' 
#' In a model predicting the first occurrence of state A, 
#' all other selected states are used as predictors. 
#' When \code{predtype} is \code{"onset-to-onset"},
#' the predictor for state B is 0 before and when that state first occurs, and 
#' it is 1 for the first state occurring directly after. 
#' Depending on \code{decaytype}, it tapers down (\code{"step"}) or 
#' remains 1 (\code{"flat"}), until a number of states have passed 
#' equal to \code{parameter}. Hence, if \code{parameter} is 3 and 
#' \code{decaytype} is "flat, it means the predictor value is 1 for 
#' the first 3 first occurrences of states after the first occurrences of state B.
#' 
#' When \code{predtype} is \code{"recent-onset-to-onset"},
#' the predictor for state B is instead 1 or tapered down from 1 to 0 
#' depending on the amount of time that has passed since 
#' the first occurrence of state B; that amount of time is given by \code{parameter}.
#' Sensible values depend on the unit of time in the date column, 
#' e.g. 24\*60\*60 for a day in seconds.
#' 
#' When \code{predtype} is \code{"recency-to-onset"}, the unit is again time,
#' but measured from the most recent occurrence of state B rather than the first.
#' Hence, if state B first occurred early on and was followed by many other states,
#' it would be 1 if it occurred again right before state A.
#' 
#' @md
#' @returns
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' 
transitionNet <- function(data,
                          predictors,
                          decaytype=c("step", "flat", "accrual"),
                          parameter,
                          decayby=c("index","date"),
                          decayfrom=c("discovery","state"),
                          direction=1,
                          sampletype="all",
                          gammas=seq(0, 1, .1),
                          alpha=1,
                          force.positive=FALSE,
                          ncores=NULL,
                          verbose=TRUE){
  predtype <- match.arg(predtype)
  mattype <- match.arg(mattype)
  decayby <- match.arg(decayby)
  decayfrom <- match.arg(decayfrom)
  
  if(decayfrom=="state"){
    if(verbose){ message("Generating predictor matrix precursor") }
    recencydata <- get_recencymat(data=data,
                                  preds=predictors,
                                  by=decayby,
                                  direction=direction,
                                  seq.onset=TRUE,
                                  verbose=verbose)
    pat <- recencydata$pat
    stairmat <- recencydata$stairmat
    rm(recencydata)
    gc()
  }
  
  if(decayfrom=="discovery"){
    stopifnot(!is.null(predictors))
    if(verbose){ message("Generating predictor matrix precursor") }
    pat <- data[!duplicated(paste(data$sequence, data$state, sep="+")), 
                c("sequence", "state", "date")]
    stairmat <- onsets2stairmat(pat=pat,
                                preds=predictors,
                                by=decayby,
                                direction=direction,
                                seq.onset=TRUE,
                                verbose=verbose)
  }
  
  # Form predictor matrix
  if(verbose){ message("Generating predictor matrix") }
  predmat <- stairmat2predmat(x=parameter, stairmat=stairmat, type=decaytype)
  rm(stairmat, data)
  gc()
  
  # Generate id matrix
  stopifnot(!is.null(predictors))
  if(verbose){ message("Generating DV matrix") }
  idmat <- onsets2idmat(pat=pat, preds=predictors)
  
  # add idx column to pat
  pat$idx <- ave(seq_along(pat$sequence),pat$sequence,FUN=function(x){seq_along(x)})
  
  if(force.positive){
    minweights <- ifelse(colnames(predmat) %in% c("(Start)","(End)","(Intercept)"),
                         -Inf, 0)
  }else{
    minweights <- rep(-Inf, ncol(predmat))
  }
  
  # Prepare cluster
  if(verbose){ message("Preparing cluster") }
  if(is.null(ncores)){ ncores <- detectCores() }
  clust <- makeCluster(ncores)
  registerDoParallel(clust, cores=ncores)
  gc(reset=T, full=T)
  on.exit({
    stopCluster(clust)
    registerDoSEQ()
  })
  
  # Run regressions
  if(verbose){ message("Running regressions on ",ncores," cores") }
  fits <- 
    foreach(currstate=colnames(idmat),
            .packages=c("glmnet","dplyr","Matrix"),
            .export=c("getStateMask","extract_coefs","EBIC.glmnet"),
            .multicombine=T) %dopar% {
              gc(reset=T, full=T)
              mask<-getStateMask(currstate=get("currstate"), pat=pat,
                                 sampletype=sampletype, direction=direction)
              
              iterfit<-glmnet(x=predmat[mask, colnames(predmat) != get("currstate")],
                              y=idmat[mask, get("currstate"), drop=F],
                              family="binomial",
                              alpha=alpha,
                              lower.limits=minweights[colnames(predmat) != get("currstate")])
              extract_coefs(iterfit, gammas)
            }
  if(verbose){ message("Assembling results") }
  names(fits) <- colnames(idmat)
  extrapreds <- ifelse(direction==1, "(Start)", "(End)")
  out <- assembleCoefficients(fits, extrapreds=extrapreds)
  return(out)
}

##########################
# Co-occurrence networks #
##########################

# This is a more efficient clone of IsingFit
cooccurrenceNet <-
  function(coocmat,
           gammas=seq(0, 1, .1),
           ncores=NULL,
           dfmax=Inf){
  coocmat <- as.matrix(coocmat)
  stopifnot(is.numeric(coocmat) & all(unique(as.vector(coocmat)) %in% 0:1))
    
  # Prepare cluster
  if(is.null(ncores)){  ncores <- detectCores() }
  clust <- makeCluster(ncores)
  registerDoParallel(clust, cores=ncores)
  on.exit({
    stopCluster(clust)
    registerDoSEQ()
  })
  
  # Run regressions
  fits <- foreach(currstate=colnames(coocmat),
                  .packages=c("glmnet", "dplyr"),
                  .export=c("extract_coefs", "EBIC.glmnet"),
                  .multicombine=T) %dopar% {
                     
    iterfit <- glmnet(x=coocmat[,colnames(coocmat) != get("currstate"), drop=F],
                      y=coocmat[,get("currstate"), drop=F],
                      family="binomial",
                      alpha=1,
                      dfmax=dfmax)
    extract_coefs(iterfit, gammas)
  }
  names(fits) <- colnames(coocmat)
  out <- assembleCoefficients(fits)
  return(out)
}
