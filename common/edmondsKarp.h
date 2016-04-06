#ifndef EDMONDS_KARP_HEADER_GUARD
#define EDMONDS_KARP_HEADER_GUARD
#include "context.h"
#include <boost/graph/filtered_graph.hpp>
#include <boost/property_map/compose_property_map.hpp>
namespace multistateTurnip
{
	struct edmondsKarpPredicate
	{
	public:
		edmondsKarpPredicate()
			:residual(NULL), graph(NULL)
		{}
		edmondsKarpPredicate(double* residual, const Context::internalDirectedGraph* graph)
			:residual(residual), graph(graph)
		{}
		typedef typename Context::internalDirectedGraph::edge_descriptor edge_descriptor;
		bool operator()(const edge_descriptor& edge) const
		{
			int index = boost::get(boost::edge_index, *graph, edge);
			return residual[index] > 0;
		}
		double* residual;
		const Context::internalDirectedGraph* graph;
	};
	struct edmondsKarpMaxFlowScratch
	{
		typedef boost::filtered_graph<Context::internalDirectedGraph, edmondsKarpPredicate> filteredGraphType;

		typedef typename boost::property_map<Context::internalDirectedGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
		typedef typename boost::property_map<Context::internalDirectedGraph, boost::edge_index_t>::const_type edgeIndexMapType;
		typedef boost::compose_property_map<double*, edgeIndexMapType> flowMapType;
		typedef boost::compose_property_map<const double*, edgeIndexMapType> capacityMapType;

		typedef typename boost::property_map<filteredGraphType, boost::vertex_index_t>::const_type filteredVertexIndexMapType;
		typedef typename boost::property_map<filteredGraphType, boost::edge_index_t>::type filteredEdgeIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<typename filteredGraphType::vertex_descriptor>::iterator, filteredVertexIndexMapType> filteredVertexPredecessorMapType;
		typedef boost::iterator_property_map<typename std::vector<int>::iterator, filteredVertexIndexMapType> filteredDistanceMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredEdgeIndexMapType> filteredColorMapType;
		typedef boost::compose_property_map<double*, filteredEdgeIndexMapType> filteredEdgeMapType;


		std::vector<filteredGraphType::vertex_descriptor> vertexPredecessor;
		std::vector<int> distance;
		std::vector<boost::default_color_type> color;
	};
	void edmondsKarpMaxFlow(const double* capacity, double* flow, double* residual, const Context::internalDirectedGraph& graph, int source, int sink, double upperBound, edmondsKarpMaxFlowScratch& scratch, double& maxFlow);
}
#endif

