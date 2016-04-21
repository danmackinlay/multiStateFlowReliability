#include "crudeMC.h"
#include "edmondsKarp.h"
namespace multistateTurnip
{
	void crudeMC(crudeMCArgs& args)
	{
		const Context& context = args.context;
		int source = context.getSource(), sink = context.getSink();

		std::vector<double>& capacityVector = context.getCapacityVector();
		
		const Context::internalGraph& undirectedGraph = context.getGraph();
		const Context::internalDirectedGraph& directedGraph = context.getDirectedGraph();
		
		std::size_t nEdges = boost::num_edges(undirectedGraph), nDirectedEdges = boost::num_edges(directedGraph);
		
		std::vector<double> flowVector(nDirectedEdges), residualVector(nDirectedEdges);
		args.count = 0;
		edmondsKarpMaxFlowScratch scratch;
		for(int i = 0; i < args.n; i++)
		{
			for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
			{
				const capacityDistribution& randomCapacityDistribution = context.getDistribution(edgeCounter);
				capacityVector[2*edgeCounter] = capacityVector[2*edgeCounter + 1] = randomCapacityDistribution(args.randomSource);
			}
			double flow = 0;
			std::fill(flowVector.begin(), flowVector.end(), 0);
			std::copy(capacityVector.begin(), capacityVector.end(), residualVector.begin());
			edmondsKarpMaxFlow(&capacityVector.front(), &flowVector.front(), &residualVector.front(), directedGraph, source, sink, args.threshold, scratch, flow);
			if(flow < args.threshold) args.count++;
		}
	}
}
