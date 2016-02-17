#include "createContext.h"
namespace multistateTurnip
{
	Context createContext(SEXP graph_sexp, capacityDistribution&& capacity, int vertex1, int vertex2, mpfr_class threshold, R_GRAPH_TYPE type)
	{
		boost::shared_ptr<std::vector<int> > interestVertices(new std::vector<int>());
		interestVertices->push_back(vertex1);
		interestVertices->push_back(vertex2);

		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph());
		convertGraph(graph_sexp, *graph, type);

		std::size_t nVertices = boost::num_vertices(*graph);
		boost::shared_ptr<std::vector<Context::vertexPosition> > vertexPositions(new std::vector<Context::vertexPosition>(nVertices, std::make_pair(0.0f, 0.0f)));
		return Context(graph, interestVertices, vertexPositions, std::move(capacity), threshold);
	}
}
