% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxSLMADS2.R
\name{coxSLMADS2}
\alias{coxSLMADS2}
\title{Fit a cox proportional hazard Model (coxph) model with pooling via Study Level Meta-Analysis (SLMA)}
\usage{
coxSLMADS2(formula, weights, dataName)
}
\arguments{
\item{formula}{a coxph() formula consistent with R syntax eg U~x+y+Z to regress
variables U on x, y and Z, where U is the survival object calculated separately and stored
in the server using the assign function ds.Surv.}

\item{weights}{a character string specifying a variable to be used as regression weights.
Specified in call to ds.coxSLMA.}

\item{dataName}{a character string specifying the name of a data.frame
holding the data for the model to be analysed under the specified model.}
}
\value{
All quantitative, Boolean, and character objects required to
enable the SLMA pooling of the separate coxph models fitted to each study -
in particular including the study-specific regression coefficients and their corresponding
standard errors.
}
\description{
This is the second serverside aggregate function called by ds.coxSLMA.
Fits a coxph model on data from single or multiple sources
with pooled co-analysis across studies being based on SLMA (Study Level Meta Analysis).
}
\details{
coxSLMADS2 is an aggregate function called by clientside function ds.coxSLMA.
ds.coxSLMA also calls another aggregate function coxSLMADS2. For more detailed
information see help for ds.coxSLMA.
}
\author{
Sofack, Ghislain.(Based on glmSLMADS2 by Paul Burton for DataSHIELD Development Team)
}
