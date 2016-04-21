#include "createContext.h"
namespace multistateTurnip
{
	Context createContext(SEXP graph_sexp, std::vector<capacityDistribution>& distributions, int vertex1, int vertex2, mpfr_class threshold, R_GRAPH_TYPE type)
	{
		std::vector<capacityDistribution> truncatedDistributions;
		for(std::vector<capacityDistribution>::iterator i = distributions.begin(); i != distributions.end(); i++)
		{
			truncatedDistributions.emplace_back(i->truncateAtMax(threshold));
		}

		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph());
		convertGraph(graph_sexp, *graph, type);

		return Context(graph, vertex1, vertex2, std::move(truncatedDistributions), threshold);
	}
}
