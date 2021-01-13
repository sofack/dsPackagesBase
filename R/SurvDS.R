#' @title Server side function called by \code{ds.Surv} to create the survival object
#' @description This function is similar to R function \code{Surv}.
#' @details This function calculates the survival object
#' usually used as a response variable in the cox model formula..
#' @param time a numeric vector indicating the start or follow up time.
#' @param time2 a numeric vector indicating the ending time of the interval for interval censored or counting process data only.
#' @param event a numeric vector providing the name of the status indicator(event).Usually binary
#' e.g 0=alive, 1=dead.For multiple endpoint data the event variable will be a factor, whose
#' first level is treated as censoring. Although unusual, the event indicator can be omitted,
#' in which case all subjects are assumed to have an event.
#' @param type character string specifying the type of censoring. Possible
#' values are "right", "left", "counting", "interval", "interval2" or "mstate"..
#' @return the results of the survival object stored on the server.
#' @author Sofack, Ghislain. (based on corTestDS by Demetris Avraam, for DataSHIELD Development Team)
#' @export
#'

SurvDS <- function (time, time2, event, type) {

  time.var <- eval(parse(text= time), envir = parent.frame())
  time2.var <- eval(parse(text=time2), envir = parent.frame())
  event.var <- eval(parse(text=event), envir = parent.frame())

  # Checking for the parameters in the parent frame

  if(!is.null(time)){
    time.var <- eval(parse(text=time), envir = parent.frame())
  }else{
    studysideMessage <- "ERROR: time must be specified as a numeric vector"
    return(list(studysideMessage=studysideMessage))
  }


  if(!is.null(time2)){
    time2.var <- eval(parse(text=time2), envir = parent.frame())
  }else{
    time2 = NULL
  }

  if(!is.null(event)){
    event.var <- eval(parse(text=event), envir = parent.frame())
  }else{
    studysideMessage <- "ERROR: event must be specified as a numeric vector"
    return(list(studysideMessage=studysideMessage))
  }

  # Evaluating the outputs based on the available inputs.

  if(!is.null(time2.var)){

    if(!is.null(type)){

      output <- survival::Surv(time = time.var, time2 = time2.var, event = event.var, type = type )
    }
    else{
      output <- survival::Surv(time = time.var, time2 = time2.var, event = event.var )
    }
  }
  else{
    if(!is.null(type)){
      output <- survival::Surv(time = time.var, event = event.var, type = type )
    }
    else{
      output <- survival::Surv(time = time.var, event = event.var )
    }
  }

}

#ASSIGN FUNCTION
# SurvDS
