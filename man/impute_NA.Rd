\name{impute_NA}
\alias{impute_NA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute Missing Values in Data Matrix}
\description{
Impute missing values in a data matrix either based on mean value or on multivariate normal distribution. 
}

\usage{
impute_NA(x, method = c("gaussian", "mean"), parameters)
}

\arguments{
  \item{x}{
Data matrix (cases in rows, variables in columns) with \code{NA} allowed.
}
  \item{method}{
'gaussian' for imputation based on multivariate normal distribution, 
'mean' for imputation based on mean vector.
}
\item{parameters}{
list of required parameters typically obtained by \code{\link{get_impute_parameters}}
}
}

\details{
For Gaussian-based imputation, this function relies on the package \code{norm} maintained by John Fox.  
}

\value{
A matrix of the same size as \code{x} with the \code{NA} replaced by their imputed values.  
}

\seealso{
\code{\link{get_impute_parameters}}
}

\examples{
## Create data with missing values 
set.seed(1234)
x <- matrix(rnorm(50), 10, 5)
x[sample(50, 10)] <- NA

## Impute missing values with mean 
pars1 <- get_impute_parameters(x, "mean")
pars1
x1 <- impute_NA(x, "mean", pars1)

## Impute missing values with Gaussian multivariate distribution
pars2 <- get_impute_parameters(x, "gaussian")
pars2
x2 <- impute_NA(x, "gaussian", pars2)
}
