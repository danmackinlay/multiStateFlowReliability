#ifndef RPACKAGE_TURNIP_HEADER_GUARD
#define RPACKAGE_TURNIP_HEADER_GUARD
#include <Rcpp.h>
namespace multistateTurnip
{
	SEXP turnip_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp, SEXP undirected_sexp);
	SEXP turnip_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp, SEXP undirected_sexp);
	SEXP turnip_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp, SEXP undirected_sexp);
}
#endif
