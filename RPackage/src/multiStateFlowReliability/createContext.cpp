#include "createContext.h"
namespace multistateTurnip
{
	Context createContext(SEXP graph_sexp, capacityDistribution& capacity, int vertex1, int vertex2, mpfr_class threshold, R_GRAPH_TYPE type)
	{
		capacityDistribution truncated = capacity.truncateAtMax(threshold);

		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph());
		convertGraph(graph_sexp, *graph, type);

		return Context(graph, vertex1, vertex2, std::move(truncated), threshold);
	}
}
