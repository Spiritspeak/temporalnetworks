
#################################
# Predictor generator functions #
#################################

onsets2idmat <- function(pat, preds){
  idmat <- vector(mode="list", length=length(preds))
  names(idmat) <- preds
  for(s in preds){
    idmat[[s]] <- Matrix(as.numeric(pat[["state"]]==s), ncol=1)
  }
  idmat <- do.call(cbind, idmat)
  colnames(idmat) <- preds
  return(idmat)
}

get_recencymat <- function(data, preds, by=c("date", "index"),
                           direction=1, seq.onset=TRUE, verbose=FALSE){
  quickave <- function(X, INDEX, FUN){
    split(X, INDEX) <- lapply(split(X, INDEX), FUN)
    X
  }
  
  data <- dplyr::arrange(data, sequence, date)
  onsets <- !duplicated(paste(data$sequence, data$state, sep="+"))
  pat <- data[onsets, c("sequence", "state", "date")]
  out <- matrix(0, ncol=length(preds), nrow=nrow(pat), dimnames=list(NULL, preds))
  idx <- 1:nrow(data)
  
  for(currstate in preds){
    if(verbose){ cat("\r", currstate, "                        ") }
    indices <- 1:nrow(data)
    indices[data$state != currstate] <- NA
    if(direction==1){
      indices <- quickave(X=indices, INDEX=data$sequence, 
                          FUN=carryforward_numeric)
    }else if(direction==-1){
      indices <- quickave(X=indices, INDEX=data$sequence, 
                          FUN=function(x){ rev(carryforward_numeric(rev(x))) })
    }
    if(by == "date"){
      lastreldate <- data$date[indices]
      out[,currstate] <- (data$date[onsets] - lastreldate[onsets]) * direction
    }else if(by == "index"){
      out[,currstate] <- (idx[onsets] - indices[onsets]) * direction
    }
  }
  rm(indices, lastreldate, onsets)
  out[is.na(out)] <- 0
  
  if(seq.onset & direction==1){
    if(verbose){ cat("\r(Start)                        ") }
    indices <- 1:nrow(data)
    sequenceonset <- lag(data$sequence) != data$sequence
    sequenceonset[1] <- T
    indices[!sequenceonset] <- NA
    indices <- carryforward_numeric(indices)
    if(by=="date"){
      out <- cbind(`(Start)`=data$date[onsets] - 
                     data$date[indices][onsets] + 1, out)
    }else if(by=="index"){
      out <- cbind(`(Start)`=idx[onsets] - 
                     indices[onsets] + 1, out)
    }
    rm(indices, sequenceonset)
  }
  
  if(seq.onset & direction==-1){
    if(verbose){ cat("\r(End)                        ") }
    indices <- 1:nrow(data)
    sequenceend <- data$sequence != lead(data$sequence)
    sequenceend[length(sequenceend)] <- T
    indices[!sequenceend] <- NA
    indices <- rev(carryforward_numeric(rev(indices)))
    if(by=="date"){
      out <- cbind(`(End)`=data$date[onsets] - 
                     data$date[indices][onsets] + 1, out)
    }else if(by=="index"){
      out <- cbind(`(End)`=idx[onsets] - 
                     indices[onsets] + 1, out)
    }
    rm(indices, sequenceend)
  }
  
  out <- Matrix(out, dimnames=list(NULL, colnames(out)))
  return(list(pat=pat, stairmat=out))
}

onsets2stairmat <- function(pat,
                            preds,
                            by=c("index", "date"),
                            direction=1,
                            seq.onset=TRUE,
                            verbose=FALSE){
  seqchains <- split(pat[["state"]], pat[["sequence"]])
  by <- match.arg(by)
  if(by == "index"){
    timechains <- split(pat[["idx"]], pat[["sequence"]])
  }else{
    timechains <- split(pat[["date"]], pat[["sequence"]])
  }
  
  allpreds <- preds
  if(seq.onset & direction == -1){ allpreds <- c("(End)", allpreds) }
  if(seq.onset & direction ==  1){ allpreds <- c("(Start)", allpreds) }
  
  stairs <- vector(length(allpreds), mode="list")
  names(stairs) <- allpreds
  for(s in preds){
    if(verbose){ cat("\r", s, "                        ") }
    currvec <- direction *
      unsplit(EnumerateFrom(sequences=seqchains,
                            times=timechains,
                            target=s),
              pat[["sequence"]])
    currvec[currvec < 0] <- 0
    stairs[[s]] <- currvec |> Matrix(ncol=1)
  }
  
  if(seq.onset & direction == 1){
    if(verbose){ cat("\r", "(Start)                        ") }
    stairs[["(Start)"]] <- 
      timechains |> 
      lapply(function(x){ x-x[1]+1 }) |> 
      unsplit(pat[["sequence"]])
  }  
  if(seq.onset & direction == -1){
    if(verbose){ cat("\r", "(End)                        ") }
    stairs[["(End)"]] <- 
      timechains |> 
      lapply(function(x){ x[length(x)]-x+1 }) |> 
      unsplit(pat[["sequence"]])
  }
  if(verbose){ cat("\r") }
  stairs <- do.call(cbind, stairs)
  colnames(stairs) <- allpreds
  return(stairs)
}


stairmat2predmat <- function(x, stairmat, 
                             type=c("step", "flat", "power", 
                                    "inverse", "accrual")){
  type <- match.arg(type)
  nonzero_stairmat <- stairmat!=0
  
  # Predictor level is x to the power of N new states entered since the state
  if(type == "power"){
    out <- Matrix(data=0,
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)))
    out[nonzero_stairmat] <- x^(stairmat[nonzero_stairmat]-1)
  }
  
  # Predictor level is 1 divided by the number of states entered since that state
  if(type == "inverse"){
    out <- Matrix(data=0,
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)))
    out[nonzero_stairmat] <- (x+1)/(x+stairmat[nonzero_stairmat])
  }
  
  # Predictor level is set to 1 for a set number of rows after a state is entered
  if(type=="flat"){
    out <- Matrix(data=as.numeric(nonzero_stairmat),
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)),
                  sparse=T)
    out[stairmat > x] <- 0
  }
  
  # Predictor level goes down in discrete steps until 0 after a state is entered
  if(type=="step"){
    nonzero_stairmat[stairmat>x] <- F
    out <- Matrix(data=0,
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)))
    out[nonzero_stairmat] <- (x+1-stairmat[nonzero_stairmat])/x
  }
  
  # Once you enter a state, the predictor is set to 1 indefinitely
  if(type=="accrual"){
    out <- nonzero_stairmat+0
  }
  
  return(out)
}

