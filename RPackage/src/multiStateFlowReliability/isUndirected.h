#ifndef IS_UNDIRECTED_HEADER_GUARD
#define IS_UNDIRECTED_HEADER_GUARD
#include <Rcpp.h>
#include "convertGraph.h"
namespace multistateTurnip
{
	bool isUndirectedIGraph(SEXP graph);
	bool isUndirectedGraphAM(SEXP graph);
	bool isUndirectedGraphNEL(SEXP graph);
	bool isUndirected(SEXP graph, R_GRAPH_TYPE graphType);
}
#endif
