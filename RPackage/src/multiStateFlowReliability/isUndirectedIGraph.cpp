#include "isUndirected.h"
#include "context.h"
namespace multistateTurnip
{
	bool isUndirectedIGraph(SEXP graph_sexp)
	{
		//Convert graph object
		Rcpp::List graph;
		try
		{
			graph = Rcpp::as<Rcpp::List>(graph_sexp);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Unable to convert input graph to a list");
		}
		bool directed = Rcpp::as<bool>(graph(1));
		return !directed;
	}
}
