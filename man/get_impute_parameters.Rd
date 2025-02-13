\name{get_impute_parameters}
\alias{get_impute_parameters}

\title{
Estimate The Parameters Of An Imputation Model For Missing Data}
\description{
Estimate the mean vector (for mean value-based imputation) or the mean vector + covariance matrix (for multivariate normal distribution-based imputation) from multivariate data with missing observations. 
}

\usage{
get_impute_parameters(x, method = c("gaussian", "mean"))
}

\arguments{
  \item{x}{
Data matrix (cases in rows, variables in columns) with \code{NA} allowed.
}
  \item{method}{
'gaussian' for imputation based on multivariate normal distribution, 
'mean' for imputation based on mean vector.
}
}

\details{
For Gaussian-based imputation, this function relies on the package \code{norm} maintained by John Fox.  
}

\value{
A list with component 'mean' if \code{method} = 'mean' ow with components 'mean' and 'sigma' (covariance matrix) if \code{method} = 'gaussian'.
}


\seealso{
\code{\link{impute_NA}}
}

\examples{
## Create data with missing values 
set.seed(1234)
x <- matrix(rnorm(50), 10, 5)
x[sample(50, 10)] <- NA

## Impute missing values with mean 
pars1 <- get_impute_parameters(x, "mean")
pars1

## Impute missing values with Gaussian multivariate distribution
pars2 <- get_impute_parameters(x, "gaussian")
pars2
}
