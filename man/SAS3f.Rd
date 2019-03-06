\name{SAS3f}
\alias{SAS3f}
\docType{data}
\title{
SAS3f database
}
\description{
A database to be used as example in the functions included on \code{unival} package. It contains the criterion, the primary factor scores and the general factor scores obtained using the program FACTOR. Those scores were obtained used a dataset of 238 responders to the Statistical Anxiety Scale (Vigil-Colet, Lorenzo-Seva, & Condon, 2008). For clarification: it does not contain the raw scores from the participant's answers to the test.
}
\usage{data("SAS3f")}
\format{
  A data frame with 238 observations and 5 variables, corresponding to the criterion, the primary factor scores and the general factor score.
}
\details{
The original test contains 24 items and measures 3 different anxiety subscales: Examination Anxiety , Asking for Help Anxiety and Interpretation Anxiety. Since they are highly correlated, they were considered related subscales from an overall scale, which measures statistical anxiety.

Since the package \code{unival} was designed for working with the factor scores and not the raw data, the provided datasets include the factor scores instead the raw data. It also contains a criterion, which in this case are the marks obtained by the responders on an Statistical exam.
}
\source{
\url{http://www.psicothema.com/PDF/3444.pdf}
}
\references{
Vigil-Colet, A., Lorenzo-Seva, U., & Condon, L. (2008). Development and validation of the Statistical Anxiety Scale. Psicothema, 20(1). \url{http://www.psicothema.com/PDF/3444.pdf}
}
\examples{
data(SAS3f)
}
\keyword{datasets}
