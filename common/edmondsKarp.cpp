#include "edmondsKarp.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>
namespace multistateTurnip
{
	double edmondsKarpMaxFlow(double* capacity, double* flow, double* residual, const Context::internalDirectedGraph& graph, int source, int sink, double upperBound, edmondsKarpMaxFlowScratch& scratch)
	{
		std::size_t nVertices = boost::num_vertices(graph);

		scratch.vertexPredecessor.resize(nVertices);
		scratch.distance.resize(nVertices);
		scratch.color.resize(nVertices);

		edmondsKarpMaxFlowScratch::edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, graph);
		edmondsKarpMaxFlowScratch::flowMapType flowMap(flow, edgeIndexMap);
		edmondsKarpMaxFlowScratch::flowMapType residualMap(residual, edgeIndexMap);
		edmondsKarpMaxFlowScratch::flowMapType capacityMap(capacity, edgeIndexMap);

		edmondsKarpPredicate predicate(residual, &graph);
		while(true)
		{
			Context::internalDirectedGraph::out_edge_iterator current, end;
			boost::tie(current, end) = boost::out_edges(source, graph);
			double currentFlow = 0;
			for(; current != end; current++)
			{
				currentFlow += flowMap[*current];
			}
			if(currentFlow > upperBound)
			{
				return currentFlow;
			}
			else
			{
				edmondsKarpMaxFlowScratch::filteredGraphType filteredGraph(graph, predicate);
				edmondsKarpMaxFlowScratch::filteredVertexIndexMapType filteredVertexIndexMap = boost::get(boost::vertex_index, filteredGraph);
				edmondsKarpMaxFlowScratch::filteredEdgeIndexMapType filteredEdgeIndexMap = boost::get(boost::edge_index, filteredGraph);
				edmondsKarpMaxFlowScratch::filteredVertexPredecessorMapType filteredVertexPredecessorMap(scratch.vertexPredecessor.begin(), filteredVertexIndexMap);
				edmondsKarpMaxFlowScratch::filteredColorMapType filteredColorMap(scratch.color.begin(), filteredEdgeIndexMap);
				edmondsKarpMaxFlowScratch::filteredDistanceMapType filteredDistanceMap(scratch.distance.begin(), filteredVertexIndexMap);
				edmondsKarpMaxFlowScratch::filteredEdgeMapType filteredResidualMap(residual, filteredEdgeIndexMap);

				boost::dijkstra_shortest_paths(filteredGraph, source, boost::predecessor_map(filteredVertexPredecessorMap).distance_map(filteredDistanceMap).weight_map(filteredResidualMap).color_map(filteredColorMap));
				//compute capacity of bottleneck
				edmondsKarpMaxFlowScratch::filteredGraphType::vertex_descriptor current = sink;
				double bottleneck = std::numeric_limits<double>::infinity();
				int bottleneckEdgeIndex = -1;
				if((int)scratch.vertexPredecessor[sink] == (int)sink)
				{
					return currentFlow;
				}
				while((int)current != source)
				{
					edmondsKarpMaxFlowScratch::filteredGraphType::vertex_descriptor next = scratch.vertexPredecessor[current];
					edmondsKarpMaxFlowScratch::filteredGraphType::edge_descriptor filteredEdge = boost::edge(next, current, filteredGraph).first;
					if(residualMap[filteredEdge] < bottleneck)
					{
						bottleneckEdgeIndex = boost::get(boost::edge_index, filteredGraph, filteredEdge);
						bottleneck = residual[bottleneckEdgeIndex];
					}
					current = next;
				}
				current = sink;
				while((int)current != source)
				{
					edmondsKarpMaxFlowScratch::filteredGraphType::vertex_descriptor next = scratch.vertexPredecessor[current];
					edmondsKarpMaxFlowScratch::filteredGraphType::edge_descriptor filteredEdge = boost::edge(next, current, filteredGraph).first;
					flowMap[filteredEdge] += bottleneck;
					residualMap[filteredEdge] -= bottleneck;

					edmondsKarpMaxFlowScratch::filteredGraphType::edge_descriptor reverseEdge = boost::get(boost::edge_reverse, filteredGraph, filteredEdge);
					flowMap[reverseEdge] = -flowMap[filteredEdge];
					residualMap[reverseEdge] = capacityMap[reverseEdge] - flowMap[reverseEdge];
					current = scratch.vertexPredecessor[current];
				}
			}
		}
	}
}
