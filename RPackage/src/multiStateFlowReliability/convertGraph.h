#ifndef CONVERT_GRAPH_HEADER_GUARD
#define CONVERT_GRAPH_HEADER_GUARD
#include <Rcpp.h>
#include "context.h"
#include "contextDirected.h"
namespace multistateTurnip
{
	enum R_GRAPH_TYPE
	{
		IGRAPH, GRAPHNEL, GRAPHAM
	};
	void convertGraph(SEXP graph_sexp, Context::inputGraph& graphRef, R_GRAPH_TYPE graphType);
	void convertGraphIGraph(SEXP graph_sexp, Context::inputGraph& graphRef);
	void convertGraphAM(SEXP graph_sexp, Context::inputGraph& graphRef);
	void convertGraphNEL(SEXP graph_sexp, Context::inputGraph& graphRef);
	
	void convertGraphDirected(SEXP graph_sexp, ContextDirected::inputGraph& graphRef, R_GRAPH_TYPE graphType);
	void convertGraphDirectedIGraph(SEXP graph_sexp, ContextDirected::inputGraph& graphRef);
	void convertGraphDirectedAM(SEXP graph_sexp, ContextDirected::inputGraph& graphRef);
	void convertGraphDirectedNEL(SEXP graph_sexp, ContextDirected::inputGraph& graphRef);
}
#endif
