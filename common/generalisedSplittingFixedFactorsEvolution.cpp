#include "generalisedSplittingFixedFactorsEvolution.h"
#include "edmondsKarp.hpp"
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "updateFlowIncremental.hpp"
namespace multistateTurnip
{
	void generalisedSplittingFixedFactorsEvolution(generalisedSplittingFixedFactorsEvolutionArgs& args)
	{
		if(args.splittingFactors.size() != args.times.size()-1)
		{
			throw std::runtime_error("Input splittingFactors must have length one less than the number of times");
		}
		if(args.times.back() != 1)
		{
			throw std::runtime_error("Last level must be 1");
		}
		args.estimate = 1;

		//Maximum number of samples at any one time.
		std::size_t productFactors = 1;
		for(std::vector<int>::iterator i = args.splittingFactors.begin(); i != args.splittingFactors.end(); i++) productFactors *= (std::size_t)*i;

		const Context& context = args.context;
		int source = context.getSource(), sink = context.getSink();
		const Context::internalGraph& graph = context.getGraph();
		const Context::internalDirectedGraph& directedGraph = context.getDirectedGraph();
		std::size_t nUndirectedEdges = boost::num_edges(graph);
		std::size_t nDirectedEdges = 2*nUndirectedEdges;

		std::vector<std::vector<double> > originalRates(nUndirectedEdges);
		std::vector<double> maximumCapacities(2*nUndirectedEdges);
		//The initial rate at the start of each PMC step
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
		std::vector<int> edgeStartingOffset(nUndirectedEdges);
		std::vector<int> nTimes(nUndirectedEdges);

		std::size_t totalTimes = 0;
		for(std::size_t i = 0; i < nUndirectedEdges; i++)
		{
			int timesThisEdge = context.getDistribution((int)i).getData().size() - 1;
			nTimes[i] = timesThisEdge;
			edgeStartingOffset[i] = totalTimes;
			totalTimes += (std::size_t)timesThisEdge;
		}

		std::vector<double> allRepairTimes(totalTimes), newAllRepairTimes;
	
		edmondsKarpMaxFlowScratch<Context::internalDirectedGraph, double> scratch;

		std::vector<double> residual((std::size_t)nDirectedEdges, 0), flow((std::size_t)nDirectedEdges, 0), capacity((std::size_t)nDirectedEdges, 0);
		std::vector<double> newResidual, newFlow, newCapacity;
		std::vector<double> maxFlows(1, 0), newMaxFlows;

		std::vector<double> tmpResidual(nDirectedEdges), tmpFlow(nDirectedEdges), tmpCapacity(nDirectedEdges);

		std::vector<int> timeCounts;
		timeCounts.resize(args.times.size(), 0);

		updateFlowIncrementalArgs<Context::internalDirectedGraph, double> updateArgs(directedGraph, scratch);
		updateArgs.nDirectedEdges = (int)nDirectedEdges;
		updateArgs.source = source;
		updateArgs.sink = sink;

		mpfr_class sum = 0, sumSquared = 0;
		for(int sampleCounter = 0; sampleCounter < args.n; sampleCounter++)
		{
			std::copy(maximumCapacities.begin(), maximumCapacities.end(), capacity.begin());
			for(int k = 0; k < (int)nUndirectedEdges; k++)
			{
				std::vector<double>& currentEdgeRates = originalRates[k];
				int timesThisEdge = nTimes[k];
				int startingOffset = edgeStartingOffset[k];
				//Simulate the anti-shock times
				for(int j = 0; j < timesThisEdge; j++)
				{
					boost::exponential_distribution<> repairDist(currentEdgeRates[j]);
					allRepairTimes[j + startingOffset] = repairDist(args.randomSource);
				}
				//Work out the actual capacaity of this edge
				const std::vector<std::pair<double, double> >& data = context.getDistribution(k).getCumulativeData();
				capacity[2*k] = capacity[2*k+1] = data.back().first;
				for(int j = 0; j < timesThisEdge; j++)
				{
					if(allRepairTimes[startingOffset + j] < args.times.front())
					{
						capacity[2*k] = capacity[2*k+1] = std::min(capacity[2*k], data[j].first);
						break;
					}
				}
			}
			std::fill(flow.begin(), flow.begin() + nDirectedEdges, 0);
			std::copy(capacity.begin(), capacity.begin()+nDirectedEdges, residual.begin());
			maxFlows[0] = 0;
			edmondsKarpMaxFlow(&(capacity[0]), &(flow[0]), &(residual[0]), directedGraph, source, sink, std::numeric_limits<double>::infinity(), scratch, maxFlows[0]);
			int currentSamples = 1;
			for(int timeCounter = 0; timeCounter < (int)args.times.size()-1; timeCounter++)
			{
				int outputCounter = 0;
				newAllRepairTimes.resize((std::size_t)args.splittingFactors[timeCounter]*(std::size_t)currentSamples * totalTimes);
				newCapacity.resize((std::size_t)args.splittingFactors[timeCounter]*(std::size_t)currentSamples * (std::size_t)nDirectedEdges);
				newResidual.resize((std::size_t)args.splittingFactors[timeCounter]*(std::size_t)currentSamples * (std::size_t)nDirectedEdges);
				newFlow.resize((std::size_t)args.splittingFactors[timeCounter]*(std::size_t)currentSamples * (std::size_t)nDirectedEdges);
				newMaxFlows.resize((std::size_t)args.splittingFactors[timeCounter]*(std::size_t)currentSamples);
				for(int sampleCounter2 = 0; sampleCounter2 < currentSamples; sampleCounter2++)
				{
					if(maxFlows[sampleCounter2] < args.level)
					{
						memcpy(&(newCapacity[outputCounter*nDirectedEdges]), &(capacity[sampleCounter2*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newFlow[outputCounter*nDirectedEdges]), &(flow[sampleCounter2*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newResidual[outputCounter*nDirectedEdges]), &(residual[sampleCounter2*nDirectedEdges]), sizeof(double)*nDirectedEdges);
						memcpy(&(newAllRepairTimes[outputCounter*totalTimes]), &(allRepairTimes[sampleCounter2*totalTimes]), sizeof(double)*totalTimes);
						newMaxFlows[outputCounter] = maxFlows[sampleCounter2];
						timeCounts[timeCounter]++;
						for(int sampleCounter3 = 0; sampleCounter3 < args.splittingFactors[timeCounter]; sampleCounter3++)
						{
							Context::internalDirectedGraph::edge_iterator edgeIterator, endEdgeIterator;
							boost::tie(edgeIterator, endEdgeIterator) = boost::edges(directedGraph);
							for(; edgeIterator != endEdgeIterator; edgeIterator++,edgeIterator++)
							{
								int edgeCounter = boost::get(boost::edge_index, directedGraph, *edgeIterator)/2;
								std::vector<double>& currentEdgeRates = originalRates[edgeCounter];
								const capacityDistribution& currentEdgeDistribution = context.getDistribution(edgeCounter);
								const std::vector<std::pair<double, double> >& currentEdgeCumulativeData = currentEdgeDistribution.getCumulativeData();
								for(int repairTimeCounter = 0; repairTimeCounter < nTimes[edgeCounter]; repairTimeCounter++)
								{
									bool unconditional = false;
									if(newAllRepairTimes[outputCounter*totalTimes + edgeStartingOffset[edgeCounter] + repairTimeCounter] >= args.times[timeCounter]) unconditional = true;
									else
									{
										for(int otherRepairTimeCounter = 0; otherRepairTimeCounter < repairTimeCounter; otherRepairTimeCounter++)
										{
											if(newAllRepairTimes[outputCounter*totalTimes + edgeStartingOffset[edgeCounter] + otherRepairTimeCounter] < newAllRepairTimes[outputCounter*totalTimes + edgeStartingOffset[edgeCounter] + repairTimeCounter])
											{
												unconditional = true; break;
											}
										}
									}
									//At this point we still have no reason to think that this part of the gibbs step should involve simulation from the unconditional distribution. We know that our resampling of this repair time might increase the time at which the max flow reaches the threshold, but does it increase by enough to make a difference?
									if(!unconditional)
									{
										std::copy(newCapacity.begin() + outputCounter*nDirectedEdges, newCapacity.begin() + (outputCounter+1)*nDirectedEdges, tmpCapacity.begin());
										std::copy(newResidual.begin() + outputCounter*nDirectedEdges, newResidual.begin() + (outputCounter+1)*nDirectedEdges, tmpResidual.begin());
										std::copy(newFlow.begin() + outputCounter*nDirectedEdges, newFlow.begin() + (outputCounter+1)*nDirectedEdges, tmpFlow.begin());
										//The new capacity of this edge if the repair time of the current anti-shock is bigger than args.times[timeCounter+1]
										double increaseTo = currentEdgeCumulativeData.back().first;
										for(int otherRepairTimeCounter = repairTimeCounter + 1; otherRepairTimeCounter < nTimes[edgeCounter]; otherRepairTimeCounter++)
										{
											if(newAllRepairTimes[outputCounter * totalTimes + edgeStartingOffset[edgeCounter] + otherRepairTimeCounter] <= args.times[timeCounter])
											{
												increaseTo = currentEdgeCumulativeData[otherRepairTimeCounter].first; break;
											}
										}
										tmpCapacity[2*edgeCounter] = tmpCapacity[2*edgeCounter + 1] = increaseTo;
										tmpResidual[2*edgeCounter] = increaseTo - tmpFlow[2*edgeCounter];
										tmpResidual[2*edgeCounter + 1] = increaseTo - tmpFlow[2*edgeCounter + 1];
										double newMaxFlow = newMaxFlows[outputCounter];
										edmondsKarpMaxFlow(&(tmpCapacity[0]), &(tmpFlow[0]), &(tmpResidual[0]), directedGraph, source, sink, std::numeric_limits<double>::infinity(), scratch, newMaxFlow);
										if(newMaxFlow < args.level) unconditional = true;
									}
									if(unconditional)
									{
										boost::exponential_distribution<> repairDist(currentEdgeRates[repairTimeCounter]);
										newAllRepairTimes[outputCounter*totalTimes + edgeStartingOffset[edgeCounter] + repairTimeCounter] = repairDist(args.randomSource);
									}
									else
									{
										boost::random::uniform_real_distribution<> uniformDist(0, 1 - exp(-currentEdgeRates[repairTimeCounter] * args.times[timeCounter]));
										newAllRepairTimes[outputCounter*totalTimes + edgeStartingOffset[edgeCounter] + repairTimeCounter] = - log(1 - uniformDist(args.randomSource)) / currentEdgeRates[repairTimeCounter];
									}
									double newCapacityThisEdge = currentEdgeCumulativeData.back().first;
									for(int otherRepairTimeCounter = 0; otherRepairTimeCounter < nTimes[edgeCounter]; otherRepairTimeCounter++)
									{
										if(newAllRepairTimes[outputCounter*totalTimes+ edgeStartingOffset[edgeCounter] + otherRepairTimeCounter] <= args.times[timeCounter])
										{
											newCapacityThisEdge = currentEdgeCumulativeData[otherRepairTimeCounter].first;
											break;
										}
									}
									updateArgs.edge = *edgeIterator;
									updateArgs.flow = &(newFlow[outputCounter*nDirectedEdges]);
									updateArgs.capacity = &(newCapacity[outputCounter*nDirectedEdges]);
									updateArgs.residual = &(newResidual[outputCounter*nDirectedEdges]);
									updateArgs.newCapacity = newCapacityThisEdge;
									double previousMaxFlow = newMaxFlows[outputCounter];
									updateFlowIncremental(updateArgs, previousMaxFlow, newMaxFlows[outputCounter]);
									if(newMaxFlows[outputCounter] < 0)
									{
										throw std::runtime_error("Internal error");
									}
								}
							}
							outputCounter++;
							if(sampleCounter3 != args.splittingFactors[timeCounter] - 1)
							{
								memcpy(&(newCapacity[outputCounter*nDirectedEdges]), &(newCapacity[(outputCounter-1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
								memcpy(&(newFlow[outputCounter*nDirectedEdges]), &(newFlow[(outputCounter-1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
								memcpy(&(newResidual[outputCounter*nDirectedEdges]), &(newResidual[(outputCounter-1)*nDirectedEdges]), sizeof(double)*nDirectedEdges);
								memcpy(&(newAllRepairTimes[outputCounter*totalTimes]), &(newAllRepairTimes[(outputCounter-1)*totalTimes]), sizeof(double)*totalTimes);
								newMaxFlows[outputCounter] = newMaxFlows[outputCounter - 1];

							}
						}
					}
				}
				capacity.swap(newCapacity);
				residual.swap(newResidual);
				flow.swap(newFlow);
				maxFlows.swap(newMaxFlows);
				allRepairTimes.swap(newAllRepairTimes);

				currentSamples = outputCounter;
				if(currentSamples == 0)
				{
					break;
				}
				//Update all the capacities, by lowing the time threshold. This removes anti-shocks, so increasing the capacities. So only one max-flow application is required. 
				for(int sampleCounter = 0; sampleCounter < currentSamples; sampleCounter++)
				{
					for(int edgeCounter = 0; edgeCounter < (int)nUndirectedEdges; edgeCounter++)
					{
						const capacityDistribution& currentEdgeDistribution = context.getDistribution(edgeCounter);
						const std::vector<std::pair<double, double> >& currentEdgeCumulativeData = currentEdgeDistribution.getCumulativeData();
						capacity[sampleCounter * nDirectedEdges + 2*edgeCounter] = capacity[sampleCounter * nDirectedEdges + 2*edgeCounter + 1] = currentEdgeCumulativeData.back().first;
						residual[sampleCounter * nDirectedEdges + 2*edgeCounter] = capacity[sampleCounter * nDirectedEdges + 2*edgeCounter] - flow[sampleCounter * nDirectedEdges + 2*edgeCounter];
						residual[sampleCounter * nDirectedEdges + 2*edgeCounter + 1] = capacity[sampleCounter * nDirectedEdges + 2*edgeCounter + 1] - flow[sampleCounter * nDirectedEdges + 2*edgeCounter + 1];
						for(int repairTimeCounter = 0; repairTimeCounter < nTimes[edgeCounter]; repairTimeCounter++)
						{
							if(allRepairTimes[sampleCounter*totalTimes + edgeStartingOffset[edgeCounter] + repairTimeCounter] < args.times[timeCounter+1])
							{
								capacity[sampleCounter * nDirectedEdges + 2*edgeCounter] = capacity[sampleCounter * nDirectedEdges + 2*edgeCounter + 1] = currentEdgeCumulativeData[repairTimeCounter].first;
								residual[sampleCounter * nDirectedEdges + 2*edgeCounter] = currentEdgeCumulativeData[repairTimeCounter].first - flow[sampleCounter * nDirectedEdges + 2*edgeCounter];
								residual[sampleCounter * nDirectedEdges + 2*edgeCounter + 1] = currentEdgeCumulativeData[repairTimeCounter].first - flow[sampleCounter * nDirectedEdges + 2*edgeCounter + 1];
								break;
							}
						}
					}
					edmondsKarpMaxFlow(&(capacity[sampleCounter*nDirectedEdges]), &(flow[sampleCounter*nDirectedEdges]), &(residual[sampleCounter*nDirectedEdges]), directedGraph, source, sink, std::numeric_limits<double>::infinity(), scratch, maxFlows[sampleCounter]);
				}
			}
			int acceptedSamples = 0;
			//Count the samples which hit the final level
			for(int sampleCounter = 0; sampleCounter < currentSamples; sampleCounter++)
			{
				if(maxFlows[sampleCounter] < args.level) 
				{
					timeCounts.back()++;
					acceptedSamples++;
				}
			}
			double ratio = (double)acceptedSamples / (double)productFactors;
			sum += ratio;
			sumSquared += ratio*ratio;
		}
		args.timeProbabilities.resize(args.times.size());
		args.timeProbabilities[0] = (double)timeCounts[0] / (double)args.n;
		for(int timeCounter = 0; timeCounter < (int)args.times.size()-1; timeCounter++)
		{
			args.timeProbabilities[timeCounter+1] = (double)timeCounts[timeCounter+1] / (double)(timeCounts[timeCounter]*args.splittingFactors[timeCounter]);
		}
		mpfr_class estimate_mpfr = sum / args.n;
		args.estimatedVariance = (double)(sumSquared/args.n - estimate_mpfr * estimate_mpfr)/args.n;
		args.estimate = (double)estimate_mpfr;
	}
}

