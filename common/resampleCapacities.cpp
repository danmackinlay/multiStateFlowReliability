#include "resampleCapacities.h"
namespace multistateTurnip
{
	void resampleCapacities(resampleCapacitiesArgs& args)
	{
		const Context::internalDirectedGraph& graph = args.context.getDirectedGraph();

		updateMaxFlowIncrementalArgs& updateMaxArgs = args.updateMaxFlowArgs;
		updateMaxArgs.nDirectedEdges = args.nDirectedEdges;
		updateMaxArgs.capacity = args.capacity;
		updateMaxArgs.flow = args.flow;
		updateMaxArgs.residual = args.residual;
		updateMaxArgs.source = args.source;
		updateMaxArgs.sink = args.sink;

		updateFlowIncrementalArgs& updateArgs = args.updateFlowArgs;
		updateArgs.source = args.source;
		updateArgs.sink = args.sink;
		updateArgs.flow = args.flow;
		updateArgs.residual = args.residual;
		updateArgs.capacity = args.capacity;
		updateArgs.nDirectedEdges = args.nDirectedEdges;

		Context::internalDirectedGraph::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(graph);
		for(; current != end; current++)
		{
			//We only want every second edge, because they come in pairs of edge / reverse edge
			current++;
			int edgeIndex = boost::get(boost::edge_index, graph, *current);
			const capacityDistribution& distribution = args.context.getDistribution(edgeIndex);
			//First work out whether there is a constraint on the flow of this edge.
			double thresholdFlowThisEdge = args.capacity[edgeIndex] + args.oldLevel - *args.maxFlow;
			double flowAfterIncrease;
			updateMaxArgs.newCapacity = thresholdFlowThisEdge;
			updateMaxArgs.edge = *current;
			updateMaxFlowIncremental(updateMaxArgs, *args.maxFlow, flowAfterIncrease);
			double newCapacity;
			//In this case we need to resample the current edge conditional on being smaller than a certain value
			if(flowAfterIncrease >= args.oldLevel)
			{
				newCapacity = distribution.sampleConditionalLessThan(args.randomSource, thresholdFlowThisEdge);
			}
			//In this case we can resample it unconditionally. 
			else
			{
				newCapacity = distribution(args.randomSource);
			}
			updateArgs.newCapacity = newCapacity;
			updateArgs.edge = *current;
			double newMaxFlow;
			updateFlowIncremental(updateArgs, *args.maxFlow, newMaxFlow);
			*args.maxFlow = newMaxFlow;
		}
	}
}
