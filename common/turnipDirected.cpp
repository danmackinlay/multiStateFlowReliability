#include "turnipDirected.h"
#include "computeConditionalProb.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <algorithm>
#include <boost/random/exponential_distribution.hpp>
#include "allPointsMaxFlow.hpp"
#include "edmondsKarp.hpp"
namespace multistateTurnip
{
	void turnipDirected(turnipDirectedArgs& args)
	{
		const ContextDirected& context = args.context;
		const ContextDirected::internalGraph& graph = context.getGraph();
		const std::size_t nVertices = boost::num_vertices(graph), nDirectedEdges = boost::num_edges(graph);
		int source = context.getSource();
		int sink = context.getSink();

		//These are only needed for the allPointsMaxFlow call. 
		allPointsMaxFlow::allPointsMaxFlowScratch<ContextDirected::internalGraph> scratch;
		//Working data for the edmonds karp call(s)
		edmondsKarpMaxFlowScratch<ContextDirected::internalGraph, double> edmondsKarpScratch;

		std::vector<double> minimumCapacities(nDirectedEdges);

		//All-points max-flow matrix
		std::vector<double> flowMatrix(nVertices*nVertices);
		//Get out the vertices for every undirected edge, in a vector. The index in the vector is the edge index of the edge
		std::vector<std::pair<ContextDirected::internalGraph::vertex_descriptor, ContextDirected::internalGraph::vertex_descriptor> > verticesPerEdge(nDirectedEdges);
		{
			ContextDirected::internalGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(graph);
			for(; current != end; current++)
			{
				int edgeIndex = boost::get(boost::edge_index, graph, *current);
				verticesPerEdge[edgeIndex] = std::make_pair(boost::source(*current, graph), boost::target(*current, graph));
			}
		}

		std::vector<std::vector<double> > originalRates(nDirectedEdges);
		std::vector<std::vector<mpfr_class> > originalRatesExact(nDirectedEdges);
		//The edges which have already been seen. This is used to exclude edges which become redundant. 
		std::vector<std::vector<bool> > alreadySeen(nDirectedEdges);
		//The rates for all the different edges, after some edges have been discarded and their rates added to some other edge
		std::vector<std::vector<mpfr_class> > ratesForEdges(nDirectedEdges);
		//The initial rate at the start of each PMC step
		mpfr_class sumAllRates = 0;
		//Set up the edge rates
		int totalLevels = 0;
		for(std::size_t i = 0; i < nDirectedEdges; i++)
		{
			std::vector<double>& currentEdgeRates = originalRates[i];
			std::vector<mpfr_class>& currentEdgeRatesExact = originalRatesExact[i];
			const capacityDistribution& currentEdgeDistribution = args.context.getDistribution(i);
			const std::vector<std::pair<double, double> >& cumulativeData = currentEdgeDistribution.getCumulativeData();
			std::size_t nLevels = currentEdgeDistribution.getData().size();
			minimumCapacities[i] = cumulativeData.front().first;
			totalLevels += nLevels - 1;
			mpfr_class cumulativeRates = 0;
			alreadySeen[i].resize(nLevels);
			ratesForEdges[i].resize(nLevels);
			for(std::size_t j = 0; j < cumulativeData.size()-1; j++)
			{
				mpfr_class newRate = -boost::multiprecision::log(mpfr_class(cumulativeData[cumulativeData.size() - j - 1].second)) - cumulativeRates;
				currentEdgeRatesExact.push_back(newRate);
				currentEdgeRates.push_back((double)newRate);
				cumulativeRates += newRate;
				sumAllRates += newRate;
			}
		}
		//The sum of all the conditional probabilities
		mpfr_class sumConditional = 0, sumSquaredConditional = 0;
		//This stores the current capacities
		std::vector<double>& capacityVector = context.getCapacityVector();
		//Residual and flow vectors for edmonds-Karp for the source and sink
		std::vector<double> flowVector(nDirectedEdges), residualVector(nDirectedEdges);
		//Residual and flow vectors for edmonds-Karp, for the max-flow across the edge which has just had its capacity increased
		std::vector<double> flowVectorIncreasedEdge(nDirectedEdges), residualVectorIncreasedEdge(nDirectedEdges);
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC;
		//Repair time vector
		std::vector<edgeRepairData> repairTimes;
		repairTimes.reserve(totalLevels);
		//Only warn about stability once
		bool warnedStability = false;

		//We store the per edge repair times seperately to start with, for the purposes of stripping out stuff that's not going to be important
		std::vector<double> perEdgeRepairTimes;

		std::vector<mpfr_class> computeConditionalProbScratch;

		//Work out the minimum possible flow
		double minimumPossibleFlow = 0;
		std::copy(minimumCapacities.begin(), minimumCapacities.end(), capacityVector.begin());
		std::copy(minimumCapacities.begin(), minimumCapacities.end(), residualVector.begin());
		std::fill(flowVector.begin(), flowVector.end(), 0);
		edmondsKarpMaxFlow<ContextDirected::internalGraph, double>(&capacityVector.front(), &flowVector.front(), &residualVector.front(), graph, source, sink, args.threshold, edmondsKarpScratch, minimumPossibleFlow);
		if(minimumPossibleFlow >= args.threshold)
		{
			args.firstMomentSingleSample = 1;
			args.secondMomentSingleSample = 1;
			args.varianceSingleSample = 0;
			args.sqrtVarianceOfEstimate = 0;
			args.relativeErrorEstimate = 0;
			return;
		}
		for(int i = 0; i < args.n; i++)
		{
			repairTimes.clear();
			//Simulate permutation via the repair times
			for(int k = 0; k < (int)nDirectedEdges; k++)
			{
				std::fill(ratesForEdges[k].begin(), ratesForEdges[k].end(), 0);

				std::vector<double>& currentEdgeRates = originalRates[k];
				std::vector<mpfr_class>& currentEdgeRatesExact = originalRatesExact[k];
				const capacityDistribution& currentEdgeDistribution = args.context.getDistribution(k);
				std::size_t nLevels = currentEdgeDistribution.getData().size();
				perEdgeRepairTimes.resize(nLevels-1);
				for(int j = 0; j < (int)nLevels - 1; j++)
				{
					boost::exponential_distribution<> repairDist(currentEdgeRates[j]);
					perEdgeRepairTimes[j] = repairDist(args.randomSource);
				}
				//The increase to highest capacity definitely occurs at some point
				edgeRepairData highest;
				highest.time = perEdgeRepairTimes[0];
				highest.rate = currentEdgeRatesExact[0];
				highest.level = 0;
				highest.edge = k;
				repairTimes.push_back(highest);

				//We can store pointers because there is no reallocation of this vector, due to reserve call
				edgeRepairData* minRepairTime = &*repairTimes.rbegin();
				for(int j = 1; j < (int)nLevels - 1; j++)
				{
					if(perEdgeRepairTimes[j] > minRepairTime->time)
					{
						minRepairTime->rate += currentEdgeRatesExact[j];
					}
					else
					{
						ratesForEdges[k][minRepairTime->level] = minRepairTime->rate;
						edgeRepairData time;
						time.time = perEdgeRepairTimes[j];
						time.level = j;
						time.edge = k;
						time.rate = currentEdgeRatesExact[j];
						repairTimes.push_back(time);
						minRepairTime = &*repairTimes.rbegin();
					}
				}
				ratesForEdges[k][minRepairTime->level] = minRepairTime->rate;
				std::fill(alreadySeen[k].begin(), alreadySeen[k].end(), false);
			}
			std::sort(repairTimes.begin(), repairTimes.end(), timeSorter);
			//No edges have yet been seen
			//The first rate is going to be sumAllRates
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			std::vector<edgeRepairData>::iterator repairTimeIterator = repairTimes.begin();
			//have we reached the point where we've got sufficient flow?
			bool insufficientFlow = true;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//The capacities and residuals are initially at the minimum possible capacitiies. The residuals are initially zero. 
			std::copy(minimumCapacities.begin(), minimumCapacities.end(), capacityVector.begin());
			std::copy(minimumCapacities.begin(), minimumCapacities.end(), residualVector.begin());
			std::fill(flowVector.begin(), flowVector.end(), 0);
			double currentFlow = 0;
			//Counter used to make sure we only call the all-points max flow once for every fixed number of steps
			int allPointsMaxFlowCounter = 1;
			while(insufficientFlow && repairTimeIterator != repairTimes.end())
			{
				//get out the parallel edge index
				int level = repairTimeIterator->level;
				//Which original edge does this correspond to?
				int edge = repairTimeIterator->edge;
				if(!alreadySeen[edge][level])
				{
					const capacityDistribution& currentEdgeDistribution = args.context.getDistribution(edge);
					const std::vector<std::pair<double, double> >& cumulativeData = currentEdgeDistribution.getCumulativeData();

					//Increase the capacity. Because we always remove edges with lower capacity, we don't need a maximum here
					double newCapacity = (cumulativeData.rbegin() + level)->first;
					double increase = newCapacity - capacityVector[edge];
					capacityVector[edge] = newCapacity;
					residualVector[edge] += increase;

					//determine whether or not we've hit the critical threshold
					edmondsKarpMaxFlow<ContextDirected::internalGraph, double>(&capacityVector.front(), &flowVector.front(), &residualVector.front(), graph, source, sink, args.threshold, edmondsKarpScratch, currentFlow);
					insufficientFlow = args.threshold > currentFlow;
					//Add the current rate
					ratesForPMC.push_back(currentRate);
					currentRate -= ratesForEdges[edge][level];
					alreadySeen[edge][level] = true;
					//If we now have newThreshold or higher flow between the the vertices for the edge,
					//we can discard ALL not yet added edges between those two vertices
					if (args.useAllPointsMaxFlow)
					{
						if ((allPointsMaxFlowCounter++ % args.allPointsMaxFlowIncrement) == 0)
						{
							allPointsMaxFlow::allPointsMaxFlow<ContextDirected::internalGraph, double>(flowMatrix, capacityVector, graph, scratch);
							ContextDirected::internalGraph::edge_iterator current, end;
							boost::tie(current, end) = boost::edges(graph);
							for (; current != end; current++)
							{
								ContextDirected::internalGraph::vertex_descriptor source = boost::source(*current, graph), target = boost::target(*current, graph);
								if (flowMatrix[source + nVertices * target] >= args.threshold)
								{
									int edgeToIgnore = boost::get(boost::edge_index, graph, *current);
									const capacityDistribution& edgeToIgnoreDistribution = args.context.getDistribution(edgeToIgnore);
									int edgeToIgnoreLevels = edgeToIgnoreDistribution.getData().size();
									for(int counter = 0; counter < edgeToIgnoreLevels; counter++)
									{
										if (!alreadySeen[edgeToIgnore][counter])
										{
											currentRate -= ratesForEdges[edgeToIgnore][counter];
											alreadySeen[edgeToIgnore][counter] = true;
										}
									}
								}
							}
						}
					}
					else 
					{
						ContextDirected::internalGraph::vertex_descriptor firstVertex = verticesPerEdge[edge].first;
						ContextDirected::internalGraph::vertex_descriptor secondVertex = verticesPerEdge[edge].second;
						//We want to start this computation from the beginning (rather than the incremental version). So the flows start as zero
						std::fill(flowVectorIncreasedEdge.begin(), flowVectorIncreasedEdge.end(), 0);
						std::copy(capacityVector.begin(), capacityVector.end(), residualVectorIncreasedEdge.begin());
						double flowIncreasedEdge = 0;
						edmondsKarpMaxFlow<ContextDirected::internalGraph, double>(&capacityVector.front(), &flowVectorIncreasedEdge.front(), &residualVectorIncreasedEdge.front(), graph, firstVertex, secondVertex, args.threshold, edmondsKarpScratch, flowIncreasedEdge);
						if(flowIncreasedEdge >= args.threshold)
						{
							const capacityDistribution& edgeToIgnoreDistribution = args.context.getDistribution(edge);
							int edgeToIgnoreLevels = edgeToIgnoreDistribution.getData().size();

							for(int counter = 0; counter < edgeToIgnoreLevels; counter++)
							{
								if(!alreadySeen[edge][counter])
								{
									currentRate -= ratesForEdges[edge][counter];
									alreadySeen[edge][counter] = true;
								}
							}
						}
					}
				}
				repairTimeIterator++;
			}
			mpfr_class additionalPart;
			//compute conditional probability
			if(!insufficientFlow)
			{
				additionalPart = computeConditionalProb(ratesForPMC, computeConditionalProbScratch);
				if(additionalPart > 1 && !warnedStability)
				{
					std::cout << "Numerical stability problem detected" << std::endl;
					warnedStability = true;
				}
			}
			else additionalPart = 1;
			//Add this conditional probability to the running sum and running sum of squares
			sumConditional += additionalPart;
			sumSquaredConditional += additionalPart*additionalPart;
		}
		args.firstMomentSingleSample = sumConditional/args.n;
		args.secondMomentSingleSample = sumSquaredConditional/args.n;
		args.varianceSingleSample = args.secondMomentSingleSample - args.firstMomentSingleSample*args.firstMomentSingleSample;
		args.sqrtVarianceOfEstimate= boost::multiprecision::sqrt(args.varianceSingleSample/args.n);
		args.relativeErrorEstimate = args.sqrtVarianceOfEstimate / args.firstMomentSingleSample;
	}
}
