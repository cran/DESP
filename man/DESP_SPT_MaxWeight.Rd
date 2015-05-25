\name{DESP_SPT_MaxWeight}
\alias{DESP_SPT_MaxWeight}
\title{Maximum weighted tree among all shortest path trees computation}
\description{
  This function selects the maximum weighted tree among all shortest path trees, for each connected components of the given graph.
}
\usage{
  DESP_SPT_MaxWeight(Graph)
}
\arguments{
  \item{Graph}{The considered graph.}
}
%\details{
%}
\value{
\code{DESP_SPT_MaxWeight} returns a list with components:
  \item{penult}{The vector containing the preceding node on the shortest path tree for each node of the graph. The root of a spanning tree is considered to be preceded by itself.}
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
