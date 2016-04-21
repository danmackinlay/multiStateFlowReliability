#ifndef CREATE_CONTEXT_HEADER_GUARD
#define CREATE_CONTEXT_HEADER_GUARD
#include "Rcpp.h"
#include "context.h"
#include "capacityDistribution.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	Context createContext(SEXP graph, std::vector<capacityDistribution>& distributions, int vertex1, int vertex2, mpfr_class threshold, R_GRAPH_TYPE graphType);
}
#endif
