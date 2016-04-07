#ifndef RESAMPLE_CAPACITIES_HEADER_GUARD
#define RESAMPLE_CAPACITIES_HEADER_GUARD
#include <vector>
#include "context.h"
#include "updateMaxFlowIncremental.h"
#include "updateFlowIncremental.h"
namespace multistateTurnip
{
	struct resampleCapacitiesArgs
	{
		resampleCapacitiesArgs(const Context& context, edmondsKarpMaxFlowScratch& scratch, boost::mt19937& randomSource)
			: context(context), scratch(scratch), updateMaxFlowArgs(context.getDirectedGraph(), scratch), updateFlowArgs(context.getDirectedGraph(), scratch), randomSource(randomSource)
		{}
		const Context& context;
		double* capacity, *flow, *residual;
		double oldLevel;
		std::vector<double> working1, working2;
		int nEdges, nDirectedEdges;
		int source, sink;
		edmondsKarpMaxFlowScratch& scratch;
		double* maxFlow;
		//In this case we just want to update the max flow
		updateMaxFlowIncrementalArgs updateMaxFlowArgs;
		//Here we want to update the max flow AND the flow and residual
		updateFlowIncrementalArgs updateFlowArgs;
		boost::mt19937& randomSource;
	};
	void resampleCapacities(resampleCapacitiesArgs& args);
}
#endif
