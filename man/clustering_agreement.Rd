\name{clustering_agreement}
\alias{clustering_agreement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measures of Agreement Between Two Clusterings
}
\description{
Return multiple measures of clustering agreement between two clusterings 
}
\usage{
clustering_agreement(x, y = NULL, method = c("chi2", "rand", "adj_rand", "fowlkes_mallow", "mirkin", "jaccard", "tpr", "fpr", "F_measure", "meila_heckerman", "max_match", "van_dongen", "mutual_info"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector representing the first clustering or a contingency table.
}
  \item{y}{
A vector of of same length as \code{x} representing the second clustering or \code{NULL} if \code{x} is a contigency table.
}
  \item{method}{
String containing the clustering measure(s) to calculate.
}
}
\details{
The clustering agreement measures may be based on concordant/discordant pairs, set overlap, or information theory. See reference for specific details. 


}
\value{
A named numeric vector containing clustering agreement measures. 
}
\references{
Wagner and Wagner (2007). Comparing Clusterings - An Overview. 
\url{https://i11www.iti.kit.edu/extra/publications/ww-cco-06.pdf}
}




\examples{
x <- sample(5, 20, replace = TRUE)
y <- sample(4, 20, replace = TRUE)
clustering_agreement(x, y)
clustering_agreement(table(x, y))
}

