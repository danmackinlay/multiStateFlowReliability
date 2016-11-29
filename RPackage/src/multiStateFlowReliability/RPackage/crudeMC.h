#ifndef RPACKAGE_CRUDEMC_HEADER_GUARD
#define RPACKAGE_CRUDEMC_HEADER_GUARD
#include <Rcpp.h>
namespace multistateTurnip
{
	SEXP crudeMC_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP undirected);
	SEXP crudeMC_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP undirected);
	SEXP crudeMC_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP undirected);
}
#endif
