#################### CREATE A BIVREC OBJECT ######################

#####
#' A function to create a bivrecSurv object
#'
#' @description
#' This function creates a bivariate recurrent response survival object (bivrecSurv object).
#'
#' @importFrom stats na.omit
#'
#' @param id Vector of subject's unique identifier (i).
#' @param episode Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param xij Vector with the lengths of time spent in event of type X for individual i in episode j.
#' @param yij vector with the lengths of time spent in event of type Y for individual i in episode j.
#' @param Ycind Vector of indicators, with values of 0 for the last episode for subject i or 1 otherwise. A subject with only one episode will have one 0.
#' @param Xcind Vector of indicators, with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have either one 1 (if he was censored at event Y) or one 0 (if he was censored at event X). A subject with censoring in event Y will have a vector of 1's.
#'
#' @return a BivRec repsonse object ready to put in a formula.
#'
#' @rdname BivRec
#' @export
#' @examples
#' set.seed(1234)
#' sim_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5))
#' bivrecsurv_data <-with(sim_data, bivrecSurv(id, epi, xij, yij, d1, d2))
#'
bivrecSurv <- function(id, episode, xij, yij, Xcind, Ycind) {

  #Check if anything is missing
  if (missing(xij)) stop("Missing - gap times for type 1 event (xij).")
  if (missing(yij)) stop("Missing - gap times for type 2 event (yij).")
  if (missing(id)) stop("Missing - subject identifiers (id).")
  if (missing(episode)) stop("Missing - episodes for each subject (episode).")
  if (missing(Ycind)) stop("Missing - censoring indicator for type 2 event (Ycind).")
  if (missing(Xcind)) stop("Missing - censoring indicator for type 1 event (Xcind).")

  #Check all vectors have same length
  all_lengths <- c(length(id),length(episode),length(xij),length(yij),length(Ycind),length(Xcind))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (id, episode, xij, yij, Ycind, Xcind) differs in length from the rest.")

  #Check xij > 0 and yij >=0 both numeric vectors
  if (!is.numeric(xij)) stop("Time arguments (xij and yij) must be numeric.")
  if (!is.numeric(yij)) stop("Time arguments (xij and yij) must be numeric.")
  if (any(xij <= 0)) stop("Time arguments for event type 1 (xij) must be positive.")
  if (any(yij < 0)) stop("Time arguments for event type 2 (yij) must be non-negative")

  #Check censoring indicators are made of only 0 or 1 values
  if (any(Xcind!=0 & Xcind!=1)) stop("Indicator vector for type 1 gap times (Xcind) must be made of 0 or 1 values only.")
  if (any(Ycind!=0 & Ycind!=1)) stop("Indicator vector for type 2 gap times (Ycind) must be made of 0 or 1 values only.")

  #ensure id's are numeric
  if (!is.numeric(id)) {
    if (is.character(id)) {id = as.numeric(as.factor(id))} else {
      if (is.factor(id)) {id = as.numeric((id))}else {
        stop("id vector must be numeric, character or factor.")}
    }
  }

  id_ref = id
  inputdf <- data.frame(id=id, epi=episode, xij=xij, yij=yij, d1=Xcind, d2=Ycind)

  #Checks for each subject
  err_xind = err_yind = err_epi = NULL
  unique_id <- unique(inputdf$id)

  for (i in 1:length(unique_id)) {
    sub_id <- unique_id[i]
    temp_by_subject <- subset(inputdf, inputdf$id==sub_id)
    temp_by_subject <- temp_by_subject[order(temp_by_subject$epi),]
    sub_n <- nrow(temp_by_subject)
    last_cx <- temp_by_subject$d1[sub_n]
    last_cy <- temp_by_subject$d2[sub_n]
    other_cx <- temp_by_subject$d1[-sub_n]
    other_cy <- temp_by_subject$d2[-sub_n]

    #Check indicators (cx is all 1's or one zero at end /  cy last is 0 and all others are 1)
    if (sum(other_cx)!=(sub_n-1)) {err_xind <- c(err_xind, sub_id)}
    if (sum(other_cy)!=(sub_n-1)) {err_yind <- c(err_yind, sub_id)}
    if (last_cy!=0) {err_yind <- c(err_yind, sub_id)}
    if (last_cx==0) {if (last_cy==1) {err_yind <- c(err_yind, sub_id)}}

    #Check episodes don't have gaps
    for (j in 1:sub_n){
      if (temp_by_subject$epi[j]!=j) {
        err_epi <- c(err_epi, sub_id)}
    }
  }

  error_subjects <- unique(c(err_xind, err_yind, err_epi))
  if (length(error_subjects>0)){
    errmsg <- paste(error_subjects, collapse = ", ")
    msg <- paste("Subjects with id", errmsg,
                 "removed becuase of gaps in episodes or incorrect values for Xcind, Ycind.",
                 sep=" ")
    print(msg)
    df4mdat <- inputdf[-which(inputdf$id %in% error_subjects), ]
  } else {df4mdat <- inputdf}

  #calculate censoring time
  ci=id2=NULL
  j=1
  df4mdat$zij <- df4mdat$xij + df4mdat$yij
  for (i in unique(df4mdat$id)){
    tempi=df4mdat[df4mdat$id == i,]
    if (nrow(tempi) == 1){
      ci=c(ci,tempi$zij)
      id2=c(id2,j)
    } else {
      ci=c(ci,rep(sum(tempi$zij),nrow(tempi)))
      id2=c(id2,rep(j,nrow(tempi)))
    }
    j=j+1
  }

  df4mdat <- cbind(id=id2, df4mdat[-1], ci)

  result <- list()
  result$id_ref = id_ref
  result$error_ids <- error_subjects
  result$data4Lreg <- mdat(dat=df4mdat) #data for Lee regression
  result$data4Creg <- df4mdat #data for Chang regression (this is also the df that is used in bivrecPlot)

  #####ADD data for cdf and marginal of NP model
  df4np <- df4mdat
  colnames(df4np)=c("id", "epi", "vij", "wij", "d1", "d2", "x0ij", "ci")
  df4np=df4np[,c("id","vij","wij","d2","d1","epi","x0ij","ci")] #change order of columns
  forcdf1 <- np.dat(df4np, ai=1)
  forcdf2 <- np.dat(df4np, ai=2)
  marg1 <- formarginal(dat = df4np) #this is from the reformat code
  marg2 <- formarginal(dat = df4np)
  formarg1 <- np.dat(marg1, ai=1)
  formarg2 <- np.dat(marg2, ai=2)
  #two np objects that have data for cdf and marg depending on ai
  result$dat4np1 <- list(forcdf=forcdf1, formarg=formarg1,refdata = df4np) #for ai=1
  result$dat4np2 <- list(forcdf=forcdf2, formarg=formarg2,refdata = df4np) #for ai=2

  class(result) <- "bivrecSurv"
  return(result)

}
