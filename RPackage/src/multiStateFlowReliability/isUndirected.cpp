#include "isUndirected.h"
namespace multistateTurnip
{
	bool isUndirected(SEXP graph_sexp, R_GRAPH_TYPE graphType)
	{
		if(graphType == IGRAPH)
		{
			return isUndirectedIGraph(graph_sexp);
		}
		else if(graphType == GRAPHAM)
		{
			return isUndirectedGraphAM(graph_sexp);
		}
		else if(graphType == GRAPHNEL)
		{
			return isUndirectedGraphNEL(graph_sexp);
		}
		else
		{
			throw std::runtime_error("Internal error");
		}
	}
}
