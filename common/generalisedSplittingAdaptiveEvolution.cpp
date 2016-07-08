#include "generalisedSplittingAdaptiveEvolution.h"
#include "edmondsKarp.hpp"
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "updateFlowIncremental.hpp"
namespace multistateTurnip
{
	namespace generalisedSplittingPrivate
	{
		struct edgeFailureData
		{
			double time;
			int edge, level;
		};
		struct networkFailureTime
		{
			double time;
			int index;
		};
		bool edgeTimeSorter(const edgeFailureData& first, const edgeFailureData& second)
		{
			return first.time < second.time;
		}
		bool networkTimeSorter(const networkFailureTime& first, const networkFailureTime& second)
		{
			return first.time < second.time;
		}
	}
	void generalisedSplittingAdaptiveEvolution(generalisedSplittingAdaptiveEvolutionArgs& args)
	{
		if(args.n % args.fraction != 0)
		{
			throw std::runtime_error("Input fraction must divide input n");
		}

		const Context& context = args.context;
		int source = context.getSource(), sink = context.getSink();
		const Context::internalGraph& graph = context.getGraph();
		const Context::internalDirectedGraph& directedGraph = context.getDirectedGraph();
		std::size_t nUndirectedEdges = boost::num_edges(graph);
		std::size_t nDirectedEdges = 2*nUndirectedEdges;

		//vector for turning edge indices into edges. 
		std::vector<Context::internalDirectedGraph::edge_descriptor> edges(nDirectedEdges);
		{
			Context::internalDirectedGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(directedGraph);
			for(; current != end; current++)
			{
				edges[boost::get(boost::edge_index, directedGraph, *current)] = *current;
			}
		}

		//The rates for each edge
		std::vector<std::vector<double> > originalRates(nUndirectedEdges);
		std::vector<double> maximumCapacities(2*nUndirectedEdges);
		for(std::size_t i = 0; i < nUndirectedEdges; i++)
		{
			std::vector<double>& currentEdgeRates = originalRates[i];
			const capacityDistribution& currentEdgeDistribution = context.getDistribution(i);
			const std::vector<std::pair<double, double> >& cumulativeData = currentEdgeDistribution.getCumulativeData();
			std::size_t nLevels = currentEdgeDistribution.getData().size();
			maximumCapacities[2*i] = maximumCapacities[2*i+1] = cumulativeData.back().first;
			mpfr_class cumulativeRates = 0;
			for(std::size_t j = 0; j < nLevels - 1; j++)
			{
				mpfr_class newRate = -boost::multiprecision::log(1 - mpfr_class(cumulativeData[j+1].second)) - cumulativeRates;
				currentEdgeRates.push_back((double)newRate);
				cumulativeRates += newRate;
			}
		}
		//Vector that gives the starting offset of the repair times for each edge
		std::vector<int> edgeStartingOffset(nUndirectedEdges);
		//Vector that gives the number of repair times, for each edge
		std::vector<int> nTimes(nUndirectedEdges);

		//Total number of repair times
		int totalTimes = 0;
		for(std::size_t i = 0; i < nUndirectedEdges; i++)
		{
			int timesThisEdge = context.getDistribution((int)i).getData().size() - 1;
			nTimes[i] = timesThisEdge;
			edgeStartingOffset[i] = totalTimes;
			totalTimes += timesThisEdge;
		}
		using namespace generalisedSplittingPrivate;

		//Repair times for each sample
		std::vector<edgeFailureData> allFailureTimes(totalTimes*args.n), newAllFailureTimes(totalTimes*args.n);
	
		edmondsKarpMaxFlowScratch<Context::internalDirectedGraph, double> scratch;

		//Max flow data at the time of failure
		std::vector<double> residual(nDirectedEdges * args.n, 0), flow(nDirectedEdges * args.n, 0), capacity(nDirectedEdges * args.n, 0);
		std::vector<double> newResidual(nDirectedEdges * args.n, 0), newFlow(nDirectedEdges * args.n, 0), newCapacity(nDirectedEdges * args.n, 0);
		std::vector<double> maxFlows(args.n, 0), newMaxFlows(args.n, 0);
		//Max flow data at the threshold time
		std::vector<double> capacityThreshold(nDirectedEdges * args.n, 0), newCapacityThreshold(nDirectedEdges * args.n, 0);
		std::vector<networkFailureTime> networkFailureTimes(args.n), newNetworkFailureTimes(args.n);
		//A vector that gives the index of each failure time within networkFailureTimes / newNetworkFailureTimes. Needed because networkFailureTimes / newNetworkFailureTimes are sorted. 
		std::vector<int> unsortedAllFailureTimes(args.n*totalTimes), newUnsortedAllFailureTimes(args.n*totalTimes);
		//Working temporaries
		std::vector<double> tmpResidual(nDirectedEdges), tmpFlow(nDirectedEdges), tmpCapacity(nDirectedEdges);

		//Work out the maximum possible max flow. 
		double maximumMaxFlow = 0;
		std::fill(tmpFlow.begin(), tmpFlow.end(), 0);
		std::vector<double> maximumFlows(nDirectedEdges), maximumResiduals(nDirectedEdges);
		std::copy(maximumCapacities.begin(), maximumCapacities.end(), maximumResiduals.begin());
		edmondsKarpMaxFlow(&(maximumCapacities[0]), &(maximumFlows[0]), &(maximumResiduals[0]), directedGraph, source, sink, std::numeric_limits<double>::infinity(), scratch, maximumMaxFlow);
		if(maximumMaxFlow < args.level) throw std::runtime_error("Internal error");

		updateFlowIncrementalArgs<Context::internalDirectedGraph, double> updateArgs(directedGraph, scratch);
		updateArgs.nDirectedEdges = (int)nDirectedEdges;
		updateArgs.source = source;
		updateArgs.sink = sink;

		//Simulate the initial samples
		for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
		{
			//Simulate all the repair times
			for(int k = 0; k < (int)nUndirectedEdges; k++)
			{
				std::vector<double>& currentEdgeRates = originalRates[k];
				int timesThisEdge = nTimes[k];
				int startingOffset = edgeStartingOffset[k];
				//Simulate the anti-shock times
				for(int j = 0; j < timesThisEdge; j++)
				{
					boost::exponential_distribution<> repairDist(currentEdgeRates[j]);
					allFailureTimes[sampleCounter * totalTimes + j + startingOffset].time = repairDist(args.randomSource);
					allFailureTimes[sampleCounter * totalTimes + j + startingOffset].level = j;
					allFailureTimes[sampleCounter * totalTimes + j + startingOffset].edge = k;
				}
			}
			//Initialise everything with the maximum possible values. 
			std::copy(maximumFlows.begin(), maximumFlows.end(), flow.begin() + sampleCounter * nDirectedEdges);
			std::copy(maximumResiduals.begin(), maximumResiduals.end(), residual.begin() + sampleCounter * nDirectedEdges);
			std::copy(maximumCapacities.begin(), maximumCapacities.end(), capacity.begin() + sampleCounter * nDirectedEdges);
			//Sort the repair times
			std::sort(allFailureTimes.begin() + sampleCounter * totalTimes, allFailureTimes.begin() + (sampleCounter + 1) * totalTimes, edgeTimeSorter);
			maxFlows[sampleCounter] = maximumMaxFlow;
			networkFailureTimes[sampleCounter].time = std::numeric_limits<double>::infinity();
			networkFailureTimes[sampleCounter].index = sampleCounter;
			for(int i = 0; i < totalTimes; i++)
			{
				edgeFailureData& relevantRepair = allFailureTimes[sampleCounter * totalTimes + i];
				int startingOffset = edgeStartingOffset[relevantRepair.edge];
				const capacityDistribution& currentEdgeDistribution = context.getDistribution(relevantRepair.edge);
				const std::vector<std::pair<double, double> >& cumulativeData = currentEdgeDistribution.getCumulativeData();

				unsortedAllFailureTimes[sampleCounter * totalTimes + startingOffset + relevantRepair.level] = i;

				if(maxFlows[sampleCounter] >= args.level)
				{
					updateArgs.capacity = &(capacity[sampleCounter*nDirectedEdges]);
					updateArgs.flow = &(flow[sampleCounter*nDirectedEdges]);
					updateArgs.residual = &(residual[sampleCounter*nDirectedEdges]);
					updateArgs.newCapacity = std::min(capacity[sampleCounter*nDirectedEdges + 2*relevantRepair.edge], cumulativeData[relevantRepair.level].first);
					updateArgs.edge = edges[2*relevantRepair.edge];
					double previousMaxFlow = maxFlows[sampleCounter];
					updateFlowIncremental(updateArgs, previousMaxFlow, maxFlows[sampleCounter]);
					networkFailureTimes[sampleCounter].time = relevantRepair.time;
				}
			}
			if(networkFailureTimes[sampleCounter].time == std::numeric_limits<double>::infinity())
			{
				throw std::runtime_error("Internal error");
			}
		}
		while(true)
		{
			std::sort(networkFailureTimes.begin(), networkFailureTimes.end(), networkTimeSorter);
			double threshold = (networkFailureTimes[args.n/args.fraction].time + networkFailureTimes[(args.n/args.fraction)-1].time)/2;
			if(threshold <= 1)
			{
				args.times.push_back(1);
				return;
			}
			args.times.push_back(threshold);
			std::fill(capacityThreshold.begin(), capacityThreshold.end(), 0);
			for(int sampleCounter = 0; sampleCounter < args.n/args.fraction; sampleCounter++)
			{
				int index = networkFailureTimes[sampleCounter].index;
				for(int i = 0; i < totalTimes; i++)
				{
					edgeFailureData& currentFailure = allFailureTimes[index*totalTimes + i];
					const std::vector<std::pair<double, double> >& cumulativeData = context.getDistribution(currentFailure.edge).getCumulativeData();
					if(currentFailure.time < threshold)
					{
						capacityThreshold[index*nDirectedEdges + 2*currentFailure.edge] = capacityThreshold[index*nDirectedEdges + 2*currentFailure.edge + 1] = std::max(capacityThreshold[index*nDirectedEdges + 2*currentFailure.edge], cumulativeData[currentFailure.level].first);
					}
				}
			}
			int outputCounter = 0;
			for(int i = 0; i < args.n/args.fraction; i++)
			{
				int index = networkFailureTimes[i].index;
				memcpy(&(newCapacity[outputCounter*nDirectedEdges]), &(capacity[index*nDirectedEdges]), sizeof(double)*nDirectedEdges);
				memcpy(&(newCapacityThreshold[outputCounter*nDirectedEdges]), &(capacityThreshold[index*nDirectedEdges]), sizeof(double)*nDirectedEdges);
				memcpy(&(newFlow[outputCounter*nDirectedEdges]), &(flow[index*nDirectedEdges]), sizeof(double)*nDirectedEdges);
				memcpy(&(newResidual[outputCounter*nDirectedEdges]), &(residual[index*nDirectedEdges]), sizeof(double)*nDirectedEdges);
				memcpy(&(newAllFailureTimes[outputCounter*totalTimes]), &(allFailureTimes[index*totalTimes]), sizeof(edgeFailureData)*totalTimes);
				memcpy(&(newUnsortedAllFailureTimes[outputCounter*totalTimes]), &(unsortedAllFailureTimes[index*totalTimes]), sizeof(int)*totalTimes);
				newMaxFlows[outputCounter] = maxFlows[index];
				newNetworkFailureTimes[outputCounter] = networkFailureTimes[i];
				for(int j = 0; j < args.fraction; j++)
				{
					networkFailureTime& currentNetworkFailure = newNetworkFailureTimes[outputCounter];
					currentNetworkFailure.index = outputCounter;

					Context::internalDirectedGraph::edge_iterator edgeIterator, endEdgeIterator;
					boost::tie(edgeIterator, endEdgeIterator) = boost::edges(directedGraph);
					for(; edgeIterator != endEdgeIterator; edgeIterator++,edgeIterator++)
					{
						int edgeCounter = boost::get(boost::edge_index, directedGraph, *edgeIterator)/2;
						std::vector<double>& currentEdgeRates = originalRates[edgeCounter];
						const capacityDistribution& currentEdgeDistribution = context.getDistribution(edgeCounter);
						const std::vector<std::pair<double, double> >& currentEdgeCumulativeData = currentEdgeDistribution.getCumulativeData();
						int startingOffset = edgeStartingOffset[edgeCounter];
						for(int repairTimeCounter = 0; repairTimeCounter < nTimes[edgeCounter]; repairTimeCounter++)
						{
							bool unconditional = false;
							int indexThisFailureTime = newUnsortedAllFailureTimes[outputCounter*totalTimes + startingOffset + repairTimeCounter];
							edgeFailureData& relevantFailure = newAllFailureTimes[outputCounter * totalTimes + indexThisFailureTime];
							double increaseTo = currentEdgeCumulativeData.back().first;
							for(int otherRepairTimeCounter = 0; otherRepairTimeCounter < nTimes[edgeCounter]; otherRepairTimeCounter++)
							{
								int indexOtherFailureTime = newUnsortedAllFailureTimes[outputCounter*totalTimes + startingOffset + otherRepairTimeCounter];
								if(newAllFailureTimes[outputCounter * totalTimes + indexOtherFailureTime].time < threshold && otherRepairTimeCounter != repairTimeCounter)
								{
									increaseTo = currentEdgeCumulativeData[otherRepairTimeCounter].first; break;
								}
							}
							if(increaseTo <= newCapacityThreshold[outputCounter*nDirectedEdges + 2*edgeCounter + 1]) unconditional = true;
							else
							{
								std::copy(newCapacityThreshold.begin() + outputCounter*nDirectedEdges, newCapacityThreshold.begin() + (outputCounter+1)*nDirectedEdges, tmpCapacity.begin());
								std::copy(newCapacityThreshold.begin() + outputCounter*nDirectedEdges, newCapacityThreshold.begin() + (outputCounter+1)*nDirectedEdges, tmpResidual.begin());
								std::fill(tmpFlow.begin(), tmpFlow.end(), 0);
								tmpCapacity[2*edgeCounter] = tmpCapacity[2*edgeCounter + 1] = increaseTo;
								tmpResidual[2*edgeCounter] = increaseTo - tmpFlow[2*edgeCounter];
								tmpResidual[2*edgeCounter + 1] = increaseTo - tmpFlow[2*edgeCounter + 1];
								double newMaxFlow = 0;
								edmondsKarpMaxFlow(&(tmpCapacity[0]), &(tmpFlow[0]), &(tmpResidual[0]), directedGraph, source, sink, std::numeric_limits<double>::infinity(), scratch, newMaxFlow);
								if(newMaxFlow < args.level) unconditional = true;
							}
							//Simulate the new failure time
							double newFailureTime;
							if(unconditional)
							{
								boost::exponential_distribution<> repairDist(currentEdgeRates[repairTimeCounter]);
								newFailureTime = repairDist(args.randomSource);
							}
							else
							{
								boost::random::uniform_real_distribution<> uniformDist(0, 1 - exp(-currentEdgeRates[repairTimeCounter] * threshold));
								newFailureTime = - log(1 - uniformDist(args.randomSource)) / currentEdgeRates[repairTimeCounter];
							}
							//Reorder the failure times (newAllFailureTimes and newUnsortedAllFailureTimes)
							edgeFailureData searchObject;
							searchObject.time = newFailureTime;
							std::vector<edgeFailureData>::iterator destination = std::lower_bound(newAllFailureTimes.begin() + outputCounter * totalTimes, newAllFailureTimes.begin() + (outputCounter + 1) * totalTimes, searchObject, edgeTimeSorter);
							std::vector<edgeFailureData>::iterator source = newAllFailureTimes.begin() + outputCounter * totalTimes + indexThisFailureTime;
							edgeFailureData copiedFailure = *source;
							copiedFailure.time = newFailureTime;
							if(source < destination)
							{
								std::copy(source+1, destination, source);
								*(destination-1) = copiedFailure;
							}
							else if(source > destination)
							{
								std::copy_backward(destination, source, source+1);
								*destination = copiedFailure;
								source++;
							}
							//update unsortedAllFailureTimes
							for(std::vector<edgeFailureData>::iterator k = std::min(source, destination); k != std::max(source, destination); k++)
							{
								newUnsortedAllFailureTimes[outputCounter*totalTimes + edgeStartingOffset[k->edge] + k->level] = (int)std::distance(newAllFailureTimes.begin() + outputCounter * totalTimes, k);
							}
							//Update the capacities so they still reflect the capacities at time currentNetworkFailure.time
							double newCapacityAfterResampling = maximumCapacities[2*edgeCounter], newCapacityAfterResamplingThreshold = maximumCapacities[2*edgeCounter];
							bool minimumFound = false, minimumFoundThreshold = false;;
							for(int k = 0; k < nTimes[edgeCounter]; k++)
							{
								int kk = newUnsortedAllFailureTimes[outputCounter*totalTimes + startingOffset + k];
								if(newAllFailureTimes[outputCounter*totalTimes + kk].time <= currentNetworkFailure.time && !minimumFound)
								{
									newCapacityAfterResampling = currentEdgeCumulativeData[k].first;
									minimumFound = true;
								}
								if(newAllFailureTimes[outputCounter*totalTimes + kk].time < threshold &&!minimumFoundThreshold)
								{
									newCapacityAfterResamplingThreshold = currentEdgeCumulativeData[k].first;
									minimumFoundThreshold = true;
								}
							}
							newCapacityThreshold[outputCounter * nDirectedEdges + 2*edgeCounter] = newCapacityThreshold[outputCounter * nDirectedEdges + 2*edgeCounter + 1] = newCapacityAfterResamplingThreshold;
							updateArgs.newCapacity = newCapacityAfterResampling;
							updateArgs.edge = edges[2*edgeCounter];
							updateArgs.capacity = &(newCapacity[outputCounter*nDirectedEdges]);
							updateArgs.flow = &(newFlow[outputCounter*nDirectedEdges]);
							updateArgs.residual = &(newResidual[outputCounter*nDirectedEdges]);
							double previousMaxFlow = newMaxFlows[outputCounter];
							updateFlowIncremental(updateArgs, previousMaxFlow, newMaxFlows[outputCounter]);
							//Now we need to iterate forwards or backwards until we find the right failure time. 
							searchObject.time = currentNetworkFailure.time;
							std::vector<edgeFailureData>::iterator currentFailureIterator = std::lower_bound(newAllFailureTimes.begin() + outputCounter * totalTimes, newAllFailureTimes.begin() + (outputCounter + 1) * totalTimes, searchObject, edgeTimeSorter);
							if(newMaxFlows[outputCounter] >= args.level)
							{
								//Iterate forwards in time. 
								while(currentFailureIterator != newAllFailureTimes.begin() + (outputCounter + 1) * totalTimes && newMaxFlows[outputCounter] >= args.level)
								{
									const capacityDistribution& currentFailureEdgeDistribution = context.getDistribution(currentFailureIterator->edge);
									updateArgs.newCapacity = std::min(currentFailureEdgeDistribution.getCumulativeData()[currentFailureIterator->level].first, newCapacity[outputCounter*nDirectedEdges + 2*currentFailureIterator->edge]);
									updateArgs.edge = edges[2*currentFailureIterator->edge];
									updateArgs.capacity = &(newCapacity[outputCounter*nDirectedEdges]);
									updateArgs.flow = &(newFlow[outputCounter*nDirectedEdges]);
									updateArgs.residual = &(newResidual[outputCounter*nDirectedEdges]);
									double previousMaxFlow = newMaxFlows[outputCounter];
									updateFlowIncremental(updateArgs, previousMaxFlow, newMaxFlows[outputCounter]);
									currentNetworkFailure.time = currentFailureIterator->time;
									currentFailureIterator++;
								}
							}
							else
							{
								//Iterate backwards in time. This means that once we have maxflow bigger than args.level, we need to step forwards one step again. 
								while(true)
								{
									const std::vector<std::pair<double, double> >& cumulativeData = context.getDistribution(currentFailureIterator->edge).getCumulativeData();
									currentNetworkFailure.time = currentFailureIterator->time;
									updateArgs.newCapacity = maximumCapacities[2*currentFailureIterator->edge];
									for(int otherRepairCounter = 0; otherRepairCounter < nTimes[currentFailureIterator->edge]; otherRepairCounter++)
									{
										int kk = newUnsortedAllFailureTimes[outputCounter * totalTimes + edgeStartingOffset[currentFailureIterator->edge] + otherRepairCounter];
										if(newAllFailureTimes[outputCounter * totalTimes + kk].time < currentFailureIterator->time)
										{
											updateArgs.newCapacity = cumulativeData[otherRepairCounter].first;
											break;
										}

									}
									if(updateArgs.newCapacity < newCapacity[outputCounter*nDirectedEdges + 2*currentFailureIterator->edge])
									{
										throw std::runtime_error("Internal error");
									}
									updateArgs.edge = edges[2*currentFailureIterator->edge];
									updateArgs.capacity = &(newCapacity[outputCounter*nDirectedEdges]);
									updateArgs.flow = &(newFlow[outputCounter*nDirectedEdges]);
									updateArgs.residual = &(newResidual[outputCounter*nDirectedEdges]);
									double previousMaxFlow = newMaxFlows[outputCounter];
									updateFlowIncremental(updateArgs, previousMaxFlow, newMaxFlows[outputCounter]);
									
									if(newMaxFlows[outputCounter] >= args.level)
									{
										updateArgs.newCapacity = cumulativeData[currentFailureIterator->level].first;
										previousMaxFlow = newMaxFlows[outputCounter];
										updateFlowIncremental(updateArgs, previousMaxFlow, newMaxFlows[outputCounter]);
										currentNetworkFailure.time = currentFailureIterator->time;
										break;
									}
									if(currentFailureIterator == newAllFailureTimes.begin() + outputCounter * totalTimes) throw std::runtime_error("Internal error");
									currentFailureIterator--;
								}
							}
						}
					}
					outputCounter++;
					if(j != args.fraction - 1)
					{
						memcpy(&(newCapacity[outputCounter*nDirectedEdges]), &(newCapacity[(outputCounter - 1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newCapacityThreshold[outputCounter*nDirectedEdges]), &(newCapacityThreshold[(outputCounter - 1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newFlow[outputCounter*nDirectedEdges]), &(newFlow[(outputCounter - 1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newResidual[outputCounter*nDirectedEdges]), &(newResidual[(outputCounter - 1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newAllFailureTimes[outputCounter*totalTimes]), &(newAllFailureTimes[(outputCounter - 1)*totalTimes]), sizeof(edgeFailureData)*totalTimes);
						memcpy(&(newUnsortedAllFailureTimes[outputCounter*totalTimes]), &(newUnsortedAllFailureTimes[(outputCounter - 1)*totalTimes]), sizeof(int)*totalTimes);
						newMaxFlows[outputCounter] = newMaxFlows[outputCounter - 1];
						newNetworkFailureTimes[outputCounter] = newNetworkFailureTimes[outputCounter - 1];
					}
				}
			}
			newCapacity.swap(capacity);
			newFlow.swap(flow);
			newResidual.swap(residual);
			newAllFailureTimes.swap(allFailureTimes);
			newUnsortedAllFailureTimes.swap(unsortedAllFailureTimes);
			newMaxFlows.swap(maxFlows);
			newNetworkFailureTimes.swap(networkFailureTimes);
		}
	}
}

