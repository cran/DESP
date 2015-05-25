\name{DESP_SPT_MaxDegreeRoot}
\alias{DESP_SPT_MaxDegreeRoot}
\title{SPT computation choosing root a priori as the node of maximal degree}
\description{
  This function computes the shortest path tree of each connected component choosing root a priori as the node of maximal degree.
}
\usage{
  DESP_SPT_MaxDegreeRoot(Graph)
}
\arguments{
  \item{Graph}{The undirected weighted graph.}
}
%\details{
%}
\value{
\code{DESP_SPT_MaxDegreeRoot} returns a list with components:
  \item{penult}{The vector containing the preceding node on the shortest path tree for each node of the graph. The root of a shortest path tree is considered to be preceded by itself.}
  \item{cc}{The list of the connected components of the given graph.}
}
\references{
  Dijkstra E. W. : A note on two problems in connexion with graphs, Numer. Math., 1, 269-271, 1959.
}

\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link[RBGL]{dijkstra.sp}}}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
