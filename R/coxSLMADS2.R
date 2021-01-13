#' @title Fit a cox proportional hazard Model (coxph) model with pooling via Study Level Meta-Analysis (SLMA)
#' @description This is the second serverside aggregate function called by ds.coxSLMA.
#' Fits a coxph model on data from single or multiple sources
#' with pooled co-analysis across studies being based on SLMA (Study Level Meta Analysis).
#' @details coxSLMADS2 is an aggregate function called by clientside function ds.coxSLMA.
#' ds.coxSLMA also calls another aggregate function coxSLMADS2. For more detailed
#' information see help for ds.coxSLMA.
#' @param formula a coxph() formula consistent with R syntax eg U~x+y+Z to regress
#' variables U on x, y and Z, where U is the survival object calculated separately and stored
#' in the server using the assign function ds.Surv.
#' @param weights a character string specifying a variable to be used as regression weights.
#' Specified in call to ds.coxSLMA.
#' @param dataName a character string specifying the name of a data.frame
#' holding the data for the model to be analysed under the specified model.
#' @return All quantitative, Boolean, and character objects required to
#' enable the SLMA pooling of the separate coxph models fitted to each study -
#' in particular including the study-specific regression coefficients and their corresponding
#' standard errors.
#' @author Sofack, Ghislain.(Based on glmSLMADS2 by Paul Burton for DataSHIELD Development Team)
#' @export


coxSLMADS2 <-function(formula, weights, dataName){


  #MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset<-as.numeric(thr$nfilter.subset)
  #nfilter.string<-as.numeric(thr$nfilter.string)


  errorMessage2<-"No errors"

  # Get the value of the 'data' parameter provided as character on the client side
  # Same is done for offset and weights lower down function


  if(!is.null(dataName)){
    dataDF <- eval(parse(text=dataName), envir = parent.frame())
  }else{
    dataDF<-NULL
  }

  # Rewrite formula extracting variables nested in structures like data frame or list
  # (e.g. D$A~D$B will be re-written A~B)
  # Note final product is a list of the variables in the model (yvector and covariates)
  # it is NOT a list of model terms - these are derived later

  # Convert formula into an editable character string
  formulatext <- Reduce(paste, deparse(formula))

  # First save original model formula

  originalFormula <- formulatext

  # # Convert formula string into separate variable names split by |
  formulatext <- gsub(" ", "", formulatext, fixed=TRUE)
  formulatext <- gsub("~", "|", formulatext, fixed=TRUE)
  formulatext <- gsub("+", "|", formulatext, fixed=TRUE)
  formulatext <- gsub("*", "|", formulatext, fixed=TRUE)
  formulatext <- gsub("||", "|", formulatext, fixed=TRUE)


  #Remember model.variables and then varnames INCLUDE BOTH yvect AND linear predictor components

  model.variables <- unlist(strsplit(formulatext, split="|", fixed=TRUE))

  varnames <- c()
  for(i in 1:length(model.variables)){
    elt <- unlist(strsplit(model.variables[i], split="$", fixed=TRUE))
    varnames <- append(varnames, elt[length(elt)])
  }

  varnames <- unique(varnames)

  if(!is.null(dataName)){
    for(v in 1:length(varnames)){
      varnames[v] <- paste0(dataName,"$",varnames[v])
      test.string.0 <- paste0(dataName,"$","0")
      test.string.1 <- paste0(dataName,"$","1")
      if(varnames[v]==test.string.0) varnames[v] <- "0"
      if(varnames[v]==test.string.1) varnames[v] <- "1"
    }
    cbindraw.text <- paste0("cbind(", paste(varnames, collapse=","), ")")
  }else{
    cbindraw.text <- paste0("cbind(", paste(varnames, collapse=","), ")")
  }


  # 	#Identify and use variable names to count missing

  all.data <- eval(parse(text=cbindraw.text), envir = parent.frame())
  all.data <- eval(parse(text=cbindraw.text))

  Ntotal <- dim(all.data)[1]

  nomiss.any <- stats::complete.cases(all.data)
  nomiss.any.data <- all.data[nomiss.any,]
  N.nomiss.any <- dim(nomiss.any.data)[1]

  Nvalid <- N.nomiss.any
  Nmissing <- Ntotal-Nvalid


  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(originalFormula))), env = parent.frame()) # here we need the formula as a 'call' object

  formula2use <- formula


  #sort out  weights

  if(is.null(weights))
  {
    varname.weights<-NULL
    cbindtext.weights <- paste0("weights.to.use <- NULL")
    eval(parse(text=cbindtext.weights), envir = parent.frame())
    weights.to.use <- NULL
  }else{
    varname.weights <- paste0(weights)
  }


  if(!(is.null(weights)))
  {
    cbindtext.weights <- paste0("weights.to.use <- cbind(", weights,")")
    eval(parse(text=cbindtext.weights), envir = parent.frame())
    cbindtext.weights <- paste0("cbind(", weights,")")
    weights.to.use <- eval(parse(text=cbindtext.weights), envir = parent.frame())
  }


  mg <- survival::coxph(formula2use, x=TRUE,  weights=weights.to.use, data=dataDF)

  y.vect<-mg$y
  X.mat<-mg$x
  pw.vect<-mg$prior.weights


  #Test for oversaturated Model

  dimX<-dim((X.mat))

  coxph.saturation.invalid<-0
  num.p<-as.numeric(dimX[2])
  num.N<-as.numeric(dimX[1])

  if(num.p>(nfilter.glm*num.N)){
    coxph.saturation.invalid <-1
    errorMessage.gos<-paste0("ERROR: Model is oversaturated (too many model parameters relative to sample size)",
                             "leading to a possible risk of disclosure - please simplify model. With ",
                             num.p," parameters and nfilter.glm = ",round(nfilter.glm,4)," you need ",
                             round((num.p/nfilter.glm),0)," observations")
  }



  #Check for Y vector validity

  y.invalid <- 0

  #Count number of unique non-missing values - disclosure risk only arises with two levels
  unique.values.noNA.y<-unique(y.vect[stats::complete.cases(y.vect)])

  #If two levels, check whether either level 0 < n < nfilter.tab
  if(length(unique.values.noNA.y)==2){
    tabvar<-table(y.vect)[table(y.vect)>=1]   #tabvar counts n in all categories with at least one observation
    min.category<-min(tabvar)
    if(min.category<nfilter.tab){
      y.invalid<-1
      errorMessage.y<-"ERROR: y (response) vector is binary with one category less than filter threshold for table cell size"
    }
  }

  #Check x matrix validity
  #Check no dichotomous x vectors with between 1 and filter.threshold
  #observations at either level

  Xpar.invalid<-rep(0,num.p)
  x.invalid<-0 #Any x parameter invalud

  for(pj in 1:num.p){
    unique.values.noNA<-unique((X.mat[,pj])[stats::complete.cases(X.mat[,pj])])

    if(length(unique.values.noNA)==2){
      tabvar<-table(X.mat[,pj])[table(X.mat[,pj])>=1] #tabvar counts n in all categories with at least one observation
      min.category<-min(tabvar)
      if(min.category<nfilter.tab){
        Xpar.invalid[pj]<-1
        x.invalid<-1
        errorMessage.x<-"ERROR: at least one column in X matrix is binary with one category less than filter threshold for table cell size"
      }
    }
  }


  #Check w vector validity
  w.invalid<-0

  if(!is.null(pw.vect))
  {
    w.vect<-pw.vect

    unique.values.noNA.w<-unique(w.vect[stats::complete.cases(w.vect)])

    if(length(unique.values.noNA.w)==2){
      tabvar<-table(w.vect)[table(w.vect)>=1]   #tabvar counts n in all categories with at least one observation
      min.category<-min(tabvar)
      if(min.category<nfilter.tab){
        w.invalid<-1
        errorMessage.w<-"ERROR: weights vector is binary with one category less than filter threshold for table cell size"
      }
    }
  }


  #Fit the model after confirming no disclosure risk

  disclosure.risk<-0

  if(y.invalid>0||w.invalid>0||sum(Xpar.invalid)>0||coxph.saturation.invalid>0){
    disclosure.risk<-1
  }

  if(disclosure.risk==0)
  {
    mg <- survival::coxph(formula2use, weights=weights.to.use, data=dataDF)

    Nvalid <- length(mg$residuals)
    Nmissing <- length(mg$na.action)
    Ntotal <- Nvalid+Nmissing

    outlist<-list( iter=mg$iter,
                   na.action=options("na.action"), call=summary(mg)$call, terms=summary(mg)$terms,
                   data=dataName, Ntotal=Ntotal, Nvalid=Nvalid, Nmissing=Nmissing,
                   weights=varname.weights,vcov= stats::vcov(mg),
                   method=mg$method, loglik=mg$loglik,
                   formula=mg$formula, coefficients=summary(mg)$coefficients)
  }else{
    errorMessage.d1<-"ERROR: Model failed in this source because of an enhanced risk of disclosure"
    errorMessage.d2<-"The following message(s) identify the cause of this enhanced risk"

    outlist.1<-list(errorMessage.1=errorMessage.d1)
    outlist.2<-list(errorMessage.2=errorMessage.d2)

    outlist.gos<-NULL
    if(coxph.saturation.invalid==1){
      outlist.gos<-list(errorMessage.gos=errorMessage.gos)
    }

    outlist.y<-NULL
    if(y.invalid==1){
      outlist.y<-list(errorMessage.y=errorMessage.y)
    }

    outlist.x<-NULL
    if(x.invalid==1){
      outlist.x<-list(errorMessage.x=errorMessage.x)
    }

    outlist.w<-NULL
    if(w.invalid==1){
      outlist.w<-list(errorMessage.w=errorMessage.w)
    }

    outlist<-list(outlist.1,outlist.2,outlist.gos,outlist.y,outlist.x,outlist.w)
  }
  #tidy up in parent.frame()
  eval(quote(rm(weights.to.use)), envir = parent.frame())
  return(outlist)

}

# AGGREGATE FUNCTION
# coxSLMADS2
