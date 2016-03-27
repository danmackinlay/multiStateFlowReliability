#ifndef UPDATE_MAX_FLOW_INCREMENTAL_HEADER_GUARD
#define UPDATE_MAX_FLOW_INCREMENTAL_HEADER_GUARD
#include "context.h"
#include "edmondsKarp.h"
namespace multistateTurnip
{
	struct updateMaxFlowIncrementalArgs
	{
		updateMaxFlowIncrementalArgs(const Context::internalDirectedGraph& graph, edmondsKarpMaxFlowScratch& scratch)
			:graph(graph), scratch(scratch)
		{}
		double* capacity, *flow, *residual;
		Context::internalDirectedGraph::edge_descriptor edge;
		const Context::internalDirectedGraph& graph;
		double newCapacity;
		edmondsKarpMaxFlowScratch& scratch;
		int source, sink;
		std::vector<double> working1, working2, working3;
		int nDirectedEdges;
	};
	void updateMaxFlowIncremental(updateMaxFlowIncrementalArgs& args, double previousMaxFlow, double& newMaxFlow);
}
#endif
