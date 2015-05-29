\name{DESP_SPT_MaxDegreeRoot2}
\alias{DESP_SPT_MaxDegreeRoot2}
\title{SPT - of maximum height equal to 1 - computation choosing root a priori as the node of maximal degree}
\description{
  This function computes iteratively the shortest path trees of maximum height equal to 1 and choosing the root a priori as the node of maximal degree, for each connected component.
}
\usage{
  DESP_SPT_MaxDegreeRoot2(Graph)
}
\arguments{
  \item{Graph}{The undirected weighted graph.}
}
%\details{
%}
\value{
\code{DESP_SPT_MaxDegreeRoot2} returns a list of lists with components (shortest path trees):
  \item{root}{The undirected weighted graph.}
  \item{penult}{The vector containing the preceding node on the shortest path tree for each node of the graph. The root of a tree is considered to be preceded by itself.}
  \item{nodes}{The set of nodes that actually belong to the tree.}
}
\references{
  Dijkstra, E. W. (1959): A note on two problems in connexion with graphs, \emph{Numer. Math.}, \bold{1}, 269--271.
}

\author{Arnak Dalalyan and Samuel Balmand.}

\seealso{\code{\link[RBGL]{dijkstra.sp}}}

%\examples{
%}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
