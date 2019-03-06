\name{unival}
\alias{unival}
\title{
Assessing essential unidimensionality using external validity information
}
\description{
Assess essential unidimensionality using external validity information.
}
\usage{
unival(y, FP, fg, PHI, FA_model = 'Lineal', type, SEP, SEG, relip,
  relig, percent = 90, display = TRUE)
}
\arguments{
  \item{y}{
Related external variable.
}
  \item{FP}{
Primary factor score estimates.
}
  \item{fg}{
General or second-order factor score estimates.
}
  \item{PHI}{
Inter-Factor correlation matrix.
}
  \item{FA_model}{
Which FA-model was used for calibration and scoring. Available options are: "Lineal" (by default) or "Graded".
}
  \item{type}{
Which type of factor score estimates were used in FP and fg. The two available options are: "ML" or "EAP" scores. If not specified, ML will be assumed.
}
  \item{SEP}{
Standard Errors (ML scores) or PSDs (EAP scores) for primary factor scores (only required when using graded model).
}
  \item{SEG}{
Standard Errors (ML scores) or PSDs (EAP scores) for the general factor (only required when when using graded model).
}
  \item{relip}{
A vector containing the marginal reliabilities of the primary factor scores estimates. It is optional except when the number of factors is 2. It can be obtained using the function fscores from the \code{mirt} package (Chalmers, 2012), or in other software like FACTOR (Lorenzo-Seva & Ferrando, 2013).
}
  \item{relig}{
The marginal reliability of the general factor (optional).
}
  \item{percent}{
Width of the confidence interval (by default 90 for 90\% confidence interval).
}
  \item{display}{
Determines if the output will be displayed in the console (TRUE by default).
}
}
\details{
\code{unival} is based on the procedure proposed by Ferrando & Lorenzo-Seva (2019). The authors proposed two group of procedures: A group of differential validity procedures to assess the extent to which the primary factor scores relate differentially to the external variables; and a group of incremental validity procedures to assess the extent to which the primary factor scores yield predictive validity increments with respect to the single general factor scores. Both groups of procedures are based on a second-order modelling schema for the general factor.

The factor scores have to be obtained externally, we suggest using FACTOR program (Lorenzo-Seva & Ferrando, 2013) or using the functions \code{\link[mirt]{mirt}}, \code{\link[mirt]{fscores}} and \code{\link[mirt]{summary-method}} included on the \code{mirt} package (Chalmers, 2012).
}
\value{
\item{differential_validity}{A vector containing the scaled disattenuated validity coefficients expected to be equal under Ho.}
\item{differential_CI}{The confidence intervals for the scaled coefficients above.}
\item{max_diffe}{The maximal difference between the most extreme scaled coefficient and the median of all of them.}
\item{maxdiffe_CI}{The confidence interval for the difference above.}
\item{contrast2}{Error corrected correlations between (a) the general factor scores and the external variable (single correlation) and (b) the multiple factor scores and the external variable (multiple correlation).}
\item{contrast2CI}{The confidence intervals for correlations above.}
\item{incremental_validity}{A value containing the difference between the single and multiple correlations above.}
\item{incremental_CI}{The confidence interval for the difference above.}
}
\references{
Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. Journal of Statistical Software, 48(6), 1-29. \url{https://doi.org/10.18637/jss.v048.i06}

Lorenzo-Seva, U., & Ferrando, P. J. (2013). Factor 9.2: A comprehensive program for fitting exploratory and semiconfirmatory factor analysis and IRT models. Applied Psychological Measurement, 37(6), 497-498. \url{https://doi.org/10.1177/0146621613487794}

Ferrando, P.J. & Lorenzo-Seva, U. (2019). An External Validity Approach for Assessing Essential Unidimensionality in Correlated-Factor Models. Educational and Psychological Measurement.

}
\author{
Pere Joan Ferrando

David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## perform unidimensionality analysis using an example dataset. The dataset is composed by the
## criterion and the factor scores, already computed using FACTOR. The correlation between factors
## was also obtained using this program. An alternative could be using the functions included in
## mirt package (Chalmers, 2012).

y = SAS3f[,1]
FP = as.matrix(SAS3f[,2:4])
fg = SAS3f[,5]
PHI = cbind(c(1,0.408,0.504),c(0.408,1,0.436),c(0.504,0.436,1))

unival(y = y, FP = FP, fg = fg, PHI = PHI)

}
