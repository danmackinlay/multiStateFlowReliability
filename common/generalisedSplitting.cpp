#include "generalisedSplitting.h"
#include "edmondsKarp.h"
#include "updateMaxFlowIncremental.h"
#include "updateFlowIncremental.h"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
namespace multistateTurnip
{
	struct resampleCapacitiesArgs
	{
		resampleCapacitiesArgs(const Context& context, edmondsKarpMaxFlowScratch& scratch, boost::mt19937& randomSource)
			: context(context), scratch(scratch), updateMaxFlowArgs(context.getDirectedGraph(), scratch), updateFlowArgs(context.getDirectedGraph(), scratch), randomSource(randomSource)
		{}
		const Context& context;
		double* capacity, *flow, *residual;
		double* newCapacity, *newFlow, *newResidual;
		double oldLevel;
		std::vector<double> working1, working2;
		int nEdges, nDirectedEdges;
		int source, sink;
		edmondsKarpMaxFlowScratch& scratch;
		double* maxFlow, *newMaxFlow;
		//In this case we just want to update the max flow
		updateMaxFlowIncrementalArgs updateMaxFlowArgs;
		//Here we want to update the max flow AND the flow and residual
		updateFlowIncrementalArgs updateFlowArgs;
		boost::mt19937& randomSource;
	};
	double resampleCapacities(resampleCapacitiesArgs& args)
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
		updateArgs.nDirectedEdges = args.nDirectedEdges;

		Context::internalDirectedGraph::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(graph);
		double maxFlow = *args.maxFlow;
		double newMaxFlow = *args.newMaxFlow;
		for(; current != end; current++)
		{
			int edgeIndex = boost::get(boost::edge_index, graph, *current);
			//First work out whether there is a constraint on the flow of this edge.
			double thresholdFlowThisEdge = args.flow[edgeIndex] + args.oldLevel - *args.maxFlow;
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
			updateFlowIncremental(updateArgs, maxFlow, newMaxFlow);
#ifndef NDEBUG
			args.working1.resize(args.nDirectedEdges);
			args.working2.resize(args.nDirectedEdges);
			memcpy(&(args.working2[0]), args.capacity, sizeof(double)*args.nDirectedEdges);
			std::fill(args.working1.begin(), args.working1.end(), 0);
			double debugValue;
			edmondsKarpMaxFlow(args.capacity, &(args.working1[0]), &(args.working2[0]), graph, args.source, args.sink, std::numeric_limits<double>::infinity(), args.scratch, debugValue);
			if(newMaxFlow != debugValue || memcmp(&(args.working1[0]), args.flow, sizeof(double)*args.nDirectedEdges) != 0)
			{
				throw std::runtime_error("Internal error");
			}
			maxFlow = newMaxFlow;
#endif
		}
		return newMaxFlow;
	}
	void generalisedSplitting(generalisedSplittingArgs& args)
	{
		args.estimate = 1;
		edmondsKarpMaxFlowScratch scratch;
		const Context::internalDirectedGraph& graph = args.context.getDirectedGraph();
		const capacityDistribution& distribution = args.context.getDistribution();
		const std::vector<int>& interestVertices = args.context.getInterestVertices();
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
		int source = interestVertices[0], sink = interestVertices[1];
		boost::random_number_generator<boost::mt19937> generator(args.randomSource);
		resampleCapacitiesArgs resampleArgs(args.context, scratch, args.randomSource);
		resampleArgs.nEdges = (int)nEdges/2;
		resampleArgs.nDirectedEdges = (int)nEdges;
		resampleArgs.source = source;
		resampleArgs.sink = sink;
#ifndef NDEBUG
		std::vector<double>& capacityVector = args.context.getCapacityVector();
#endif
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
			edmondsKarpMaxFlow(currentSampleCapacity, currentSampleFlow, currentSampleResidual, graph, source, sink, oldLevel, scratch, maxFlows[sampleCounter]);
#ifndef NDEBUG
			double flow2 = args.context.getMaxFlow(capacityVector);
			if(maxFlows[sampleCounter] != flow2) throw std::runtime_error("Internal error");
#endif
			if(maxFlows[sampleCounter] < oldLevel)
			{
				samplesLessThanFlow++;
			}
		}
		args.estimate = (double)samplesLessThanFlow / (double)args.n;
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
					for(int copyCounter = 0; copyCounter < copiesThisSample; copyCounter++)
					{
						resampleArgs.capacity = &(capacity[sampleCounter*nEdges]);
						resampleArgs.residual = &(residual[sampleCounter*nEdges]);
						resampleArgs.flow = &(flow[sampleCounter*nEdges]);
						resampleArgs.newCapacity = &(newCapacity[outputCounter*nEdges]);
						resampleArgs.newResidual = &(newResidual[outputCounter*nEdges]);
						resampleArgs.newFlow = &(newFlow[outputCounter*nEdges]);
						resampleArgs.oldLevel = oldLevel;
						resampleArgs.maxFlow = &(maxFlows[sampleCounter]);
						resampleArgs.newMaxFlow = &(newMaxFlows[outputCounter]);
						resampleCapacities(resampleArgs);
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

			copies = args.n / samplesLessThanFlow;
			extra = args.n - copies * samplesLessThanFlow;
			for(int i = 0; i < extra; i++) extraCopy[shuffled[i]] = true;

			//Update the estimate
			args.estimate *= (double)samplesLessThanFlow / (double)args.n;
			oldLevel = newLevel;
		}
	}
}
