#include "generalisedSplittingFixedFactors.h"
#include "edmondsKarp.h"
#include "updateMaxFlowIncremental.h"
#include "updateFlowIncremental.h"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <iostream>
#include "resampleCapacities.h"
namespace multistateTurnip
{
	void generalisedSplittingFixedFactors(generalisedSplittingFixedFactorsArgs& args)
	{
		if(args.splittingFactors.size() != args.levels.size()-1)
		{
			throw std::runtime_error("Input splittingFactors must have length one less than the number of levels");
		}
		args.estimate = 1;

		//Maximum number of samples at any one time. 
		int productFactors = 1;
		for(std::vector<int>::iterator i = args.splittingFactors.begin(); i != args.splittingFactors.end(); i++) productFactors *= *i;

		edmondsKarpMaxFlowScratch scratch;
		const Context::internalDirectedGraph& graph = args.context.getDirectedGraph();
		int source = args.context.getSource(), sink = args.context.getSink();
		std::size_t nEdges = boost::num_edges(graph);
		std::vector<double> residual(nEdges * productFactors, 0), flow(nEdges * productFactors, 0), capacity(nEdges * productFactors, 0);
		std::vector<double> newResidual(nEdges * productFactors, 0), newFlow(nEdges * productFactors, 0), newCapacity(nEdges * productFactors);
		int samplesLessThanFlow = 0;
		std::vector<double> maxFlows(productFactors, 0), newMaxFlows(productFactors, 0);

		std::vector<int> levelCounts;
		levelCounts.resize(args.levels.size(), 0);

		resampleCapacitiesArgs resampleArgs(args.context, scratch, args.randomSource);
		resampleArgs.nEdges = (int)nEdges/2;
		resampleArgs.nDirectedEdges = (int)nEdges;
		resampleArgs.source = source;
		resampleArgs.sink = sink;
		//Generate initial sample
		for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
		{
			double oldLevel = args.levels.front();
			double* currentSampleCapacity = &capacity.front();
			double* currentSampleResidual = &residual.front();
			double* currentSampleFlow = &flow.front();
			for(int edgeCounter = 0; edgeCounter < (int)(nEdges/2); edgeCounter++)
			{
				const capacityDistribution& distribution = args.context.getDistribution(edgeCounter);
				currentSampleResidual[2*edgeCounter] = currentSampleResidual[2*edgeCounter+1] = currentSampleCapacity[2*edgeCounter] = currentSampleCapacity[2*edgeCounter+1] = distribution(args.randomSource);
			}
			maxFlows[0] = 0;
			std::fill(currentSampleFlow, currentSampleFlow + nEdges, 0);
			edmondsKarpMaxFlow(currentSampleCapacity, currentSampleFlow, currentSampleResidual, graph, source, sink, oldLevel, scratch, maxFlows[0]);
			int currentSamples = 1;
			for(int levelCounter = 0; levelCounter < (int)args.levels.size()-1; levelCounter++)
			{
				double oldLevel = args.levels[levelCounter];
				int outputCounter = 0;
				for(int sampleCounter = 0; sampleCounter < currentSamples; sampleCounter++)
				{
					if(maxFlows[sampleCounter] < oldLevel)
					{
						levelCounts[levelCounter]++;
						memcpy(&(newCapacity[outputCounter*nEdges]), &(capacity[sampleCounter*nEdges]), sizeof(double)*nEdges);
						memcpy(&(newFlow[outputCounter*nEdges]), &(flow[sampleCounter*nEdges]), sizeof(double)*nEdges);
						memcpy(&(newResidual[outputCounter*nEdges]), &(residual[sampleCounter*nEdges]), sizeof(double)*nEdges);
						newMaxFlows[outputCounter] = maxFlows[sampleCounter];
						for(int sampleCounter2 = 0; sampleCounter2 < args.splittingFactors[levelCounter]; sampleCounter2++)
						{
							resampleArgs.capacity = &(newCapacity[outputCounter*nEdges]);
							resampleArgs.residual = &(newResidual[outputCounter*nEdges]);
							resampleArgs.flow = &(newFlow[outputCounter*nEdges]);
							resampleArgs.oldLevel = oldLevel;
							resampleArgs.maxFlow = &(newMaxFlows[outputCounter]);
							resampleCapacities(resampleArgs);
							if(sampleCounter2 != args.splittingFactors[levelCounter] - 1)
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

				currentSamples = outputCounter;
				if(currentSamples == 0)
				{
					break;
				}
			}
			//Count the samples which hit the final level
			for(int sampleCounter = 0; sampleCounter < currentSamples; sampleCounter++)
			{
				if(maxFlows[sampleCounter] < args.levels.back()) levelCounts.back()++;
			}
		}
		args.levelProbabilities.resize(args.levels.size());
		args.estimate = args.levelProbabilities[0] = (double)levelCounts[0] / (double)args.n;
		for(int levelCounter = 0; levelCounter < (int)args.levels.size()-1; levelCounter++)
		{
			args.levelProbabilities[levelCounter+1] = (double)levelCounts[levelCounter+1] / (double)(levelCounts[levelCounter]*args.splittingFactors[levelCounter]);
			args.estimate *= args.levelProbabilities[levelCounter+1];
		}
	}
}
