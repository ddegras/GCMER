\name{transitive_closure}
\alias{transitive_closure}

\title{
Determine The Transitive Closure Of A Set Of Edges In A Graph
}

\description{
Given a subset of edges in a graph and an optional list of nodes, this function adds all required edges such that the graph has fully connected components (= cliques). In the context of binary relations, this function finds the transitive closure of a homogeneous binary relation. 
}

\usage{
transitive_closure(edges, nodes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{
Two-column matrix of class \code{numeric} or \code{character} 
specifying a subset of edges in a graph or a homogeneous binary relation.
}
  \item{nodes}{
Optional vector containing all nodes in the graph or the set over which the relation is defined. 
}
}

\details{
If provided, \code{nodes} must be of same class as \code{edges} and must contain all the values in \code{edges}.

The binary relation is assumed to be reflexive and symmetric. As such, its transitive closure forms an equivalence relation.
}

\value{
A LIST with fields
\item{edges}{The completed 2-column matrix such that all connected components are cliques (or all pairs form a transitive closure). Loops are not included.}
\item{comp_vec}{integer vector specifying the component/equivalence class of each node.}
\item{ncomps}{number of graph components (=cliques)}
}


\seealso{
\code{\link{find_connected_components}}
}

\examples{
# Consider graph with one circular component (1,2,3) and one singleton (4)
# Only list edges 1-2 and 2-3
nodes <- 1:4
edges <- matrix(c(1,2,2,3), 2, 2, byrow = TRUE)
transitive_closure(edges, nodes)
}

