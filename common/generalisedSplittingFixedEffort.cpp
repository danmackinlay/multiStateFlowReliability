#include "generalisedSplittingFixedEffort.h"
#include "edmondsKarp.h"
#include "updateMaxFlowIncremental.h"
#include "updateFlowIncremental.h"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include <iostream>
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
	void resampleCapacities(resampleCapacitiesArgs& args)
	{
		const capacityDistribution& capacity = args.context.getDistribution();
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
				newCapacity = capacity.sampleConditionalLessThan(args.randomSource, thresholdFlowThisEdge);
			}
			//In this case we can resample it unconditionally. 
			else
			{
				newCapacity = capacity(args.randomSource);
			}
			updateArgs.newCapacity = newCapacity;
			updateArgs.edge = *current;
			double newMaxFlow;
			updateFlowIncremental(updateArgs, *args.maxFlow, newMaxFlow);
			*args.maxFlow = newMaxFlow;
		}
	}
	void generalisedSplittingFixedEffort(generalisedSplittingFixedEffortArgs& args)
	{
		args.estimate = 1;
		edmondsKarpMaxFlowScratch scratch;
		const Context::internalDirectedGraph& graph = args.context.getDirectedGraph();
		const capacityDistribution& distribution = args.context.getDistribution();
		int source = args.context.getSource(), sink = args.context.getSink();
		std::size_t nEdges = boost::num_edges(graph);
		std::vector<double> residual(nEdges * args.n, 0), flow(nEdges * args.n, 0), capacity(nEdges * args.n, 0);
		std::vector<double> newResidual(nEdges * args.n, 0), newFlow(nEdges * args.n, 0), newCapacity(nEdges * args.n, 0);
		int samplesLessThanFlow = 0;
		std::vector<double> maxFlows(args.n, 0), newMaxFlows(args.n, 0);
		//Used for the stratified sampling
		std::vector<int> shuffled;
		shuffled.reserve(args.n);
		std::vector<bool> extraCopy(args.n, false);
		double oldLevel = args.levels.front();
		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		resampleCapacitiesArgs resampleArgs(args.context, scratch, args.randomSource);
		resampleArgs.nEdges = (int)nEdges/2;
		resampleArgs.nDirectedEdges = (int)nEdges;
		resampleArgs.source = source;
		resampleArgs.sink = sink;
		//Generate initial sample
		for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
		{
			double* currentSampleCapacity = &(capacity[sampleCounter*nEdges]);
			double* currentSampleResidual = &(residual[sampleCounter*nEdges]);
			double* currentSampleFlow = &(flow[sampleCounter*nEdges]);
			for(int edgeCounter = 0; edgeCounter < (int)(nEdges/2); edgeCounter++)
			{
				currentSampleResidual[2*edgeCounter] = currentSampleResidual[2*edgeCounter+1] = currentSampleCapacity[2*edgeCounter] = currentSampleCapacity[2*edgeCounter+1] = distribution(args.randomSource);
			}
			edmondsKarpMaxFlow(currentSampleCapacity, currentSampleFlow, currentSampleResidual, graph, source, sink, oldLevel, scratch, maxFlows[sampleCounter]);
			if(maxFlows[sampleCounter] < oldLevel)
			{
				samplesLessThanFlow++;
			}
		}
		args.estimate = (double)samplesLessThanFlow / (double)args.n;
		if(samplesLessThanFlow == 0)
		{
			args.estimate = 0;
			return;
		}
		int copies = args.n / samplesLessThanFlow;
		int extra = args.n - copies * samplesLessThanFlow;
		for(std::vector<double>::iterator i = maxFlows.begin(); i != maxFlows.end(); i++)
		{
			if(*i < oldLevel) shuffled.push_back((int)std::distance(maxFlows.begin(), i));
		}
		boost::range::random_shuffle(shuffled, generator);
		for(int i = 0; i < extra; i++) extraCopy[shuffled[i]] = true;

		shuffled.clear();
		for(int levelCounter = 1; levelCounter < (int)args.levels.size(); levelCounter++)
		{
			double newLevel = args.levels[levelCounter];
			int outputCounter = 0;
			for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
			{
				if(maxFlows[sampleCounter] < oldLevel)
				{
					int copiesThisSample = copies;
					if(extraCopy[sampleCounter]) copiesThisSample++;
					memcpy(&(newCapacity[outputCounter*nEdges]), &(capacity[sampleCounter*nEdges]), sizeof(double)*nEdges);
					memcpy(&(newFlow[outputCounter*nEdges]), &(flow[sampleCounter*nEdges]), sizeof(double)*nEdges);
					memcpy(&(newResidual[outputCounter*nEdges]), &(residual[sampleCounter*nEdges]), sizeof(double)*nEdges);
					newMaxFlows[outputCounter] = maxFlows[sampleCounter];
					for(int copyCounter = 0; copyCounter < copiesThisSample; copyCounter++)
					{
						resampleArgs.capacity = &(newCapacity[outputCounter*nEdges]);
						resampleArgs.residual = &(newResidual[outputCounter*nEdges]);
						resampleArgs.flow = &(newFlow[outputCounter*nEdges]);
						resampleArgs.oldLevel = oldLevel;
						resampleArgs.maxFlow = &(newMaxFlows[outputCounter]);
						resampleCapacities(resampleArgs);
						if(copyCounter != copiesThisSample-1)
						{
							memcpy(&(newCapacity[(outputCounter+1)*nEdges]), &(newCapacity[outputCounter*nEdges]), sizeof(double)*nEdges);
							memcpy(&(newResidual[(outputCounter+1)*nEdges]), &(newResidual[outputCounter*nEdges]), sizeof(double)*nEdges);
							memcpy(&(newFlow[(outputCounter+1)*nEdges]), &(newFlow[outputCounter*nEdges]), sizeof(double)*nEdges);
							newMaxFlows[outputCounter+1] = newMaxFlows[outputCounter];
						}
						outputCounter++;
					}
				}
			}
			capacity.swap(newCapacity);
			residual.swap(newResidual);
			flow.swap(newFlow);
			maxFlows.swap(newMaxFlows);
			//Now work out the number of samples
			std::fill(extraCopy.begin(), extraCopy.end(), false);
			samplesLessThanFlow = 0;
			shuffled.clear();
			for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
			{
				if(maxFlows[sampleCounter] <  newLevel)
				{
					samplesLessThanFlow++;
					shuffled.push_back(sampleCounter);
				}
			}
			boost::range::random_shuffle(shuffled, generator);

			if(samplesLessThanFlow == 0)
			{
				args.estimate = 0; 
				return;
			}
			copies = args.n / samplesLessThanFlow;
			extra = args.n - copies * samplesLessThanFlow;
			for(int i = 0; i < extra; i++) extraCopy[shuffled[i]] = true;

			//Update the estimate
			args.estimate *= (double)samplesLessThanFlow / (double)args.n;

			oldLevel = newLevel;
		}
	}
}
