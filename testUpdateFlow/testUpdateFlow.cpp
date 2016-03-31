#include "updateFlowIncremental.h"
#include "updateMaxFlowIncremental.h"
#include <limits>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/graph/connected_components.hpp>
#include "edmondsKarp.h"
namespace multistateTurnip
{
	int main(int argc, char **argv)
	{
		edmondsKarpMaxFlowScratch scratch;

		//construct random graph
		boost::mt19937 randomSource;
		int edgeIndex = 0;
		const int nVertices = 20;

		boost::random::bernoulli_distribution<> bern;

		Context::internalDirectedGraph randomGraph(nVertices);
		for(int i = 0; i < nVertices; i++)
		{
			for(int j = 0; j < i; j++)
			{
				if(bern(randomSource))
				{
					Context::internalDirectedGraph::edge_descriptor newEdge = boost::add_edge(i, j, randomGraph).first;
					Context::internalDirectedGraph::edge_descriptor newReverseEdge = boost::add_edge(j, i, randomGraph).first;
					boost::put(boost::edge_index, randomGraph, newEdge, edgeIndex);
					edgeIndex++;
					boost::put(boost::edge_index, randomGraph, newReverseEdge, edgeIndex);
					edgeIndex++;
					boost::put(boost::edge_reverse, randomGraph, newEdge, newReverseEdge);
					boost::put(boost::edge_reverse, randomGraph, newReverseEdge, newEdge);
				}
			}
		}
		std::vector<int> componentMap;
		std::size_t nDirectedEdges = boost::num_edges(randomGraph);
		componentMap.resize(nDirectedEdges);
		int components = boost::connected_components(randomGraph, &(componentMap[0]));
		if(components != 1) throw std::runtime_error("Random graph did not have a single component");

		//The edmonds karp variables
		std::vector<double> capacity(nDirectedEdges, 0), flow(nDirectedEdges, 0), residual(nDirectedEdges, 0);
		//The incremental flow variables
		std::vector<double> incCapacity(nDirectedEdges, 0), incFlow(nDirectedEdges, 0), incResidual(nDirectedEdges, 0);
		//distribution for the edges
		boost::random::exponential_distribution<double> expDist(3);

		//incremental arguments
		edmondsKarpMaxFlowScratch incrementalScratch;

		updateFlowIncrementalArgs incrementalFlowArgs(randomGraph, incrementalScratch);
		incrementalFlowArgs.capacity = &(incCapacity[0]);
		incrementalFlowArgs.flow = &(incFlow[0]);
		incrementalFlowArgs.residual = &(incResidual[0]);
		incrementalFlowArgs.source = 0;
		incrementalFlowArgs.sink = nVertices - 1;
		incrementalFlowArgs.nDirectedEdges = nDirectedEdges;

		updateMaxFlowIncrementalArgs incrementalMaxFlowArgs(randomGraph, incrementalScratch);
		incrementalMaxFlowArgs.capacity = &(incCapacity[0]);
		incrementalMaxFlowArgs.flow = &(incFlow[0]);
		incrementalMaxFlowArgs.residual = &(incResidual[0]);
		incrementalMaxFlowArgs.source = 0;
		incrementalMaxFlowArgs.sink = nVertices - 1;
		incrementalMaxFlowArgs.nDirectedEdges = nDirectedEdges;

		double incrementalMaxFlow = 0;
		//Change all the edges 10000 times
		for(int i = 0; i < 10000; i++)
		{
			Context::internalDirectedGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(randomGraph);
			for(; current != end; current++)
			{
				int edgeIndex = boost::get(boost::edge_index, randomGraph, *current);
				Context::internalDirectedGraph::edge_descriptor reverseEdge = boost::get(boost::edge_reverse, randomGraph, *current);
				int reverseEdgeIndex = boost::get(boost::edge_index, randomGraph, reverseEdge);

				double newCapacity = expDist(randomSource);
				capacity[edgeIndex] = capacity[reverseEdgeIndex] = newCapacity;
				std::fill(flow.begin(), flow.end(), 0);
				std::copy(capacity.begin(), capacity.end(), residual.begin());
				double edmondsKarpMaxFlowResult = 0;
				edmondsKarpMaxFlow(&(capacity[0]), &(flow[0]), &(residual[0]), randomGraph, 0, nVertices - 1, std::numeric_limits<double>::infinity(), scratch, edmondsKarpMaxFlowResult);

				//Test the code that only updates the max flow
				incrementalMaxFlowArgs.edge = *current;
				incrementalMaxFlowArgs.newCapacity = newCapacity;
				double previousMaxFlow = incrementalMaxFlow;
				double onlyMaxFlow;
				updateMaxFlowIncremental(incrementalMaxFlowArgs, previousMaxFlow, onlyMaxFlow);
				if(fabs(onlyMaxFlow - edmondsKarpMaxFlowResult) > 1e-8)
				{
					throw std::runtime_error("Max flow mismatch");
				}

				//Test the code that updates everything
				incrementalFlowArgs.edge = *current;
				incrementalFlowArgs.newCapacity = newCapacity;
				updateFlowIncremental(incrementalFlowArgs, previousMaxFlow, incrementalMaxFlow);
				if(fabs(incrementalMaxFlow - edmondsKarpMaxFlowResult) > 1e-8)
				{
					throw std::runtime_error("Max flow mismatch");
				}

				//Check that the incremental code has worked. 
				double tmp = 0;
				edmondsKarpMaxFlow(&(incCapacity[0]), &(incFlow[0]), &(incResidual[0]), randomGraph, 0, nVertices - 1, std::numeric_limits<double>::infinity(), scratch, tmp);
				if(fabs(tmp) > 1e-8)
				{
					throw std::runtime_error("Max flow mismatch");
				}
				//The flows and residuals don't have to be the same. 
				/*if(memcmp(&(flow[0]), &(incFlow[0]), sizeof(double)*nDirectedEdges) != 0)
				{
					throw std::runtime_error("Flow mismatch");
				}
				if(memcmp(&(residual[0]), &(incResidual[0]), sizeof(double)*nDirectedEdges) != 0)
				{
					throw std::runtime_error("Residual mismatch");
				}*/

			}
		}
		std::cout << "Success!" << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return multistateTurnip::main(argc, argv);
}
