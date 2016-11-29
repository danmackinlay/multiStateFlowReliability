#ifndef RPACKAGE_PMC_HEADER_GUARD
#define RPACKAGE_PMC_HEADER_GUARD
#include <Rcpp.h>
namespace multistateTurnip
{
	SEXP pmc_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP undirected, SEXP verbose);
	SEXP pmc_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP undirected, SEXP verbose);
	SEXP pmc_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices, SEXP undirected, SEXP verbose);
}
#endif
