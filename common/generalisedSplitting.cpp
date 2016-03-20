#include "generalisedSplitting.h"
#include "edmondsKarp.h"
namespace multistateTurnip
{
	void generalisedSplitting(generalisedSplittingArgs& args)
	{
		edmondsKarpMaxFlowScratch scratch;
		const Context::internalDirectedGraph& graph = args.context.getDirectedGraph();
		const capacityDistribution& distribution = args.context.getDistribution();
		const std::vector<int>& interestVertices = args.context.getInterestVertices();
		std::size_t nEdges = boost::num_edges(graph);
		std::vector<double> residual(nEdges * args.n, 0), flow(nEdges * args.n, 0), capacity(nEdges * args.n, 0);
		int counter = 0;
		std::vector<double>& capacityVector = args.context.getCapacityVector();
		//Generate initial sample
		for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
		{
			double* currentSampleCapacity = &(capacity[sampleCounter*nEdges]);
			double* currentSampleResidual = &(residual[sampleCounter*nEdges]);
			double* currentSampleFlow = &(flow[sampleCounter*nEdges]);
			for(int edgeCounter = 0; edgeCounter < (int)(nEdges/2); edgeCounter++)
			{
#ifndef NDEBUG
				capacityVector[2*edgeCounter] = capacityVector[2*edgeCounter+1] = 
#endif
				currentSampleResidual[2*edgeCounter] = currentSampleResidual[2*edgeCounter+1] = currentSampleCapacity[2*edgeCounter] = currentSampleCapacity[2*edgeCounter+1] = distribution(args.randomSource);
			}
			double flow = edmondsKarpMaxFlow(currentSampleCapacity, currentSampleFlow, currentSampleResidual, graph, interestVertices[0], interestVertices[1], args.levels.back(), scratch);
#ifndef NDEBUG
			double flow2 = args.context.getMaxFlow(capacityVector);
			if(flow != flow2) throw std::runtime_error("Internal error");
#endif
			if(flow < args.levels.back()) counter++;
		}
		args.estimate = (double)counter / (double)args.n;
	}
}
