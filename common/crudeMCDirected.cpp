#include "crudeMCDirected.h"
#include "edmondsKarp.hpp"
namespace multistateTurnip
{
	void crudeMCDirected(crudeMCDirectedArgs& args)
	{
		const ContextDirected& context = args.context;
		int source = context.getSource(), sink = context.getSink();

		std::vector<double>& capacityVector = context.getCapacityVector();
		
		const ContextDirected::internalGraph& directedGraph = context.getGraph();
		
		std::size_t nDirectedEdges = boost::num_edges(directedGraph);
		
		std::vector<double> flowVector(nDirectedEdges), residualVector(nDirectedEdges);
		args.count = 0;
		edmondsKarpMaxFlowScratch<ContextDirected::internalGraph, double> scratch;
		for(int i = 0; i < args.n; i++)
		{
			for(std::size_t edgeCounter = 0; edgeCounter < nDirectedEdges; edgeCounter++)
			{
				const capacityDistribution& randomCapacityDistribution = context.getDistribution(edgeCounter);
				capacityVector[edgeCounter] = randomCapacityDistribution(args.randomSource);
			}
			double flow = 0;
			std::fill(flowVector.begin(), flowVector.end(), 0);
			std::copy(capacityVector.begin(), capacityVector.end(), residualVector.begin());
			edmondsKarpMaxFlow<ContextDirected::internalGraph, double>(&capacityVector.front(), &flowVector.front(), &residualVector.front(), directedGraph, source, sink, args.threshold, scratch, flow);
			if(flow < args.threshold) args.count++;
		}
	}
}
