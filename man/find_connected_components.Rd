\name{find_connected_components}
\alias{find_connected_components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find The Connected Components of an Adjacency Matrix
}
\description{
Find the connected components of an adjacency matrix
}
\usage{
find_connected_components(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A square matrix containing TRUE/FALSE or numerical values. 
}
}
\details{
The entries of the matrix \code{x} will be converted to logical type. The function checks if  \code{x} is symmetric and throws an error message if not.  


}
\value{
An integer vector of length \code{nrow(x)} containing the component numbers (1, 2, ...)
}

\examples{
x <- matrix(FALSE, 10, 10)
x[1:2,1:2] <- TRUE
x[3:6,3:6] <- TRUE
x[7:9,7:9] <- TRUE
find_connected_components(x)
}
