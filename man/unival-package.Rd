\name{unival-package}
\alias{unival-package}
\title{
Assessing essential unidimensionality using external validity information
}
\description{
Package for assessing the unidimensionality of a set of items using external validity information. It can be applied on linear or graded factor analytic models.
}
\details{ \code{unival} is based on the procedure proposed by Ferrando & Lorenzo-Seva (2019). The authors proposed two group of procedures: A group of differential validity procedures to assess the extent to which the primary factor scores relate differentially to the external variables; and a group of incremental validity procedures to assess the extent to which the primary factor scores yield predictive validity increments with respect to the single general factor scores. Both groups of procedures are based on a second-order modelling schema for the general factor.

The factor scores have to be obtained externally, we suggest using FACTOR program (Lorenzo-Seva & Ferrando, 2013) or using the functions \code{\link[mirt]{mirt}}, \code{\link[mirt]{fscores}} and \code{\link[mirt]{summary-method}} included on the \code{mirt} package (Chalmers, 2012).

More information can be found on the documentation page of the function \code{\link{unival}}.

}
\value{
\item{\code{\link{unival}}}{Assess essential unidimensionality using external validity information.}
}
\references{
Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory Package for the R Environment. Journal of Statistical Software, 48(6), 1-29. \doi{10.18637/jss.v048.i06}

Lorenzo-Seva, U., & Ferrando, P. J. (2013). Factor 9.2: A comprehensive program for fitting exploratory and semiconfirmatory factor analysis and IRT models. Applied Psychological Measurement, 37(6), 497-498. \doi{10.1177/0146621613487794}

Ferrando, P.J. & Lorenzo-Seva, U. (2019). An External Validity Approach for Assessing Essential Unidimensionality in Correlated-Factor Models. Educational and Psychological Measurement. \doi{10.1177/0013164418824755}
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
## \code{mirt} package (Chalmers, 2012).

y=SAS3f[,1]
FP=as.matrix(SAS3f[,2:4])
fg=SAS3f[,5]
PHI=cbind(c(1,0.408,0.504),c(0.408,1,0.436),c(0.504,0.436,1))

unival(y = y, FP = FP, fg = fg, PHI = PHI)
}
\keyword{ package }
