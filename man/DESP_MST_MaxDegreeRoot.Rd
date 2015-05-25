\name{DESP_MST_MaxDegreeRoot}
\alias{DESP_MST_MaxDegreeRoot}
\title{MST computation}
\description{
  This function computes the minimum spanning trees, choosing the root of each tree as the node of maximal degree of the connected component.
}
\usage{
  DESP_MST_MaxDegreeRoot(Graph)
}
\arguments{
  \item{Graph}{The considered graph.}
}
%\details{
%}
\value{
\code{DESP_MST_MaxDegreeRoot} returns a list with components:
  \item{penult}{The vector containing the preceding node on the spanning tree for each node of the graph. The root of a spanning tree is considered to be preceded by itself.}
  \item{cc}{The list of the connected components of the given graph.}
}

\references{
  Kruskal Joseph B., Jr. : On the shortest spanning subtree of a graph and the traveling salesman problem, Proc. Amer. Math. Soc., 7, 48-50, 1956.
}

\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link[graph]{connComp}}, \code{\link[RBGL]{mstree.kruskal}}}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
