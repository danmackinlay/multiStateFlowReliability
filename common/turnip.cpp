#include "turnip.h"
#include "computeConditionalProb.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <algorithm>
#include <boost/random/exponential_distribution.hpp>
#include "allPointsMaxFlow.hpp"
namespace multistateTurnip
{
	namespace turnipPrivate
	{
		bool secondArgumentSorter(const std::pair<int, double>& first, const std::pair<int, double>& second)
		{
			return first.second < second.second;
		}
	}
	void turnip(turnipArgs& args)
	{
		const Context& context = args.context;
		const capacityDistribution& distribution = context.getDistribution();
		const std::vector<std::pair<double, double> >& cumulativeData = distribution.getCumulativeData();
		std::size_t nLevels = distribution.getData().size();
		std::size_t nParallelEdges = context.getNEdges() * (cumulativeData.size()-1);

		//Get out the vertices for every edge, in a vector. 
		std::vector<std::pair<Context::internalGraph::vertex_descriptor, Context::internalGraph::vertex_descriptor> > verticesPerEdge;
		const Context::internalGraph& graph = context.getGraph();
		//These are only needed for the allPointsMaxFlow call. Otherwise for the single max-flow call we use Context.maxFlow, which
		//uses the directed graph internally. 
		const Context::internalDirectedGraph& directedGraph = context.getDirectedGraph();
		allPointsMaxFlow::allPointsMaxFlowScratch<Context::internalDirectedGraph> scratch;
		const std::size_t nVertices = boost::num_vertices(graph);
		std::vector<double> flowMatrix(nVertices*nVertices);
		{
			Context::internalGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(graph);
			for(; current != end; current++)
			{
				verticesPerEdge.push_back(std::make_pair(current->m_source, current->m_target));
			}
		}

		//Describe all the parrallel edges - In terms of the rate of that parallel edge, 
		//which original edge it belongs to, and which index of parallel edge it is among all the parallel edges 
		//for that original edge
		std::vector<int> originalEdgeIndex;
		std::vector<int> originalEdgeLevel;
		std::vector<double> originalRates;
		//Set up the original edge indices and rates
		for(std::size_t i = 0; i < context.getNEdges(); i++)
		{
			double cumulativeRates = 0;
			originalEdgeIndex.insert(originalEdgeIndex.end(), cumulativeData.size()-1, (int)i);
			originalEdgeLevel.insert(originalEdgeLevel.end(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)cumulativeData.size()-1));
			for(std::size_t j = 0; j < cumulativeData.size()-1; j++)
			{
				double newRate = -log(cumulativeData[cumulativeData.size() - j - 1].second) - cumulativeRates;
				originalRates.push_back(newRate);
				cumulativeRates += newRate;
			}
		}
		boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::sum_kahan> > acc;
		std::for_each(originalRates.begin(), originalRates.end(), std::ref(acc));
		//The sum of all the conditional probabilities
		mpfr_class sumConditional = 0, sumSquaredConditional = 0;
		//The initial rate at the start of each PMC step
		double sumAllRates = boost::accumulators::sum_kahan(acc);
		//This stores the current capacities
		std::vector<double>& capacityVector = context.getCapacityVector();
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC;
		//Repair time vector
		std::vector<std::pair<int, double> > repairTimes;
		repairTimes.resize(nParallelEdges);
		//The edges which have already been seen. This is used to exclude edges which become redundant. 
		std::vector<bool> alreadySeen(nParallelEdges);
		//Only warn about stability once
		bool warnedStability = false;
		for(int i = 0; i < args.n; i++)
		{
			//Simulate permutation
			for(std::size_t j = 0; j < nParallelEdges; j++)
			{
				boost::exponential_distribution<> repairDist(originalRates[j]);
				repairTimes[j].second = repairDist(args.randomSource);
				repairTimes[j].first = (int)j;
			}
			std::sort(repairTimes.begin(), repairTimes.end(), turnipPrivate::secondArgumentSorter);
			//No edges have yet been seen
			std::fill(alreadySeen.begin(), alreadySeen.end(), false);
			//The first rate is going to be this
			double currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			int permutationCounter = 0;
			//have we reached the point where we've got sufficient flow?
			bool insufficientFlow = true;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//The capacities are initially zero
			std::fill(capacityVector.begin(), capacityVector.end(), 0);
			//Counter used to make sure we only call the all-points max flow once for every fixed number of steps
			int allPointsMaxFlowCounter = 1;
			while(insufficientFlow)
			{
				//get out the parallel edge index
				int parallelEdgeIndex = repairTimes[permutationCounter].first;
				//Which original edge does this correspond to?
				int originalEdgeIndexThisLoop = originalEdgeIndex[parallelEdgeIndex];
				if(!alreadySeen[parallelEdgeIndex])
				{
					//Increase the capacity. Because we always remove edges with lower capacity, we don't need a maximum here
					capacityVector[2*originalEdgeIndexThisLoop] = capacityVector[2*originalEdgeIndexThisLoop+1] = (cumulativeData.rbegin() + originalEdgeLevel[parallelEdgeIndex])->first;
					//determine whether or not we've hit the critical threshold
					double currentFlow = context.getMaxFlow(capacityVector);
					insufficientFlow = args.threshold > currentFlow;
					//Add the current rate
					ratesForPMC.push_back(currentRate);
					
					//Start going back through the other parallel edges for this original edge
					//if they hadn't already been observed to occur, adjust the rate to indicate that we don't need them.
					do
					{
						if (originalEdgeIndex[parallelEdgeIndex] != originalEdgeIndexThisLoop) break;
						if (!alreadySeen[parallelEdgeIndex])
						{
							currentRate -= originalRates[parallelEdgeIndex];
							alreadySeen[parallelEdgeIndex] = true;
						}
						parallelEdgeIndex++;
					} while (parallelEdgeIndex >= 0 && parallelEdgeIndex < (int)nParallelEdges);

					//If we now have newThreshold or higher flow between the the vertices for edge originalEdgeIndexThisLoop,
					//we can discard ALL not yet added edges between those two vertices
					Context::internalGraph::vertex_descriptor firstVertex = verticesPerEdge[originalEdgeIndexThisLoop].first;
					Context::internalGraph::vertex_descriptor secondVertex = verticesPerEdge[originalEdgeIndexThisLoop].second;
					if (args.useAllPointsMaxFlow)
					{
						if ((allPointsMaxFlowCounter++ % args.allPointsMaxFlowIncrement) == 0)
						{
							allPointsMaxFlow::allPointsMaxFlow<Context::internalDirectedGraph>(flowMatrix, capacityVector, directedGraph, scratch);
							Context::internalGraph::edge_iterator current, end;
							boost::tie(current, end) = boost::edges(graph);
							for (; current != end; current++)
							{
								Context::internalGraph::vertex_descriptor source = current->m_source, target = current->m_target;
								if (flowMatrix[source + nVertices * target] >= args.threshold)
								{
									originalEdgeIndexThisLoop = boost::get(boost::edge_index, graph, *current);
									int parallelEdgeIndex = (int)(originalEdgeIndexThisLoop*(nLevels - 1));
									do
									{
										if (originalEdgeIndex[parallelEdgeIndex] != originalEdgeIndexThisLoop) break;
										if (!alreadySeen[parallelEdgeIndex])
										{
											currentRate -= originalRates[parallelEdgeIndex];
											alreadySeen[parallelEdgeIndex] = true;
										}
										parallelEdgeIndex++;
									} while (parallelEdgeIndex >= 0 && parallelEdgeIndex < (int)nParallelEdges);
								}
							}
						}
					}
					else if(context.getMaxFlow(capacityVector, firstVertex, secondVertex) >= args.threshold)
					{
						parallelEdgeIndex = (int)(originalEdgeIndexThisLoop*(nLevels-1));
						do
						{
							if(originalEdgeIndex[parallelEdgeIndex] != originalEdgeIndexThisLoop) break;
							if(!alreadySeen[parallelEdgeIndex]) 
							{
								currentRate -= originalRates[parallelEdgeIndex];
								alreadySeen[parallelEdgeIndex] = true;
							}
							parallelEdgeIndex++;
						}
						while(parallelEdgeIndex >= 0 && parallelEdgeIndex < (int)nParallelEdges);
					}
				}
				permutationCounter++;
			}
			mpfr_class additionalPart = computeConditionalProb(ratesForPMC);
			if(additionalPart > 1 && !warnedStability)
			{
				std::cout << "Numerical stability problem detected" << std::endl;
				warnedStability = true;
			}
			sumConditional += additionalPart;
			sumSquaredConditional += additionalPart*additionalPart;
		}
		args.estimateFirstMoment = sumConditional/args.n;
		args.estimateSecondMoment = sumSquaredConditional/args.n;
		args.varianceEstimate = args.estimateSecondMoment - args.estimateFirstMoment*args.estimateFirstMoment;
		args.sqrtVarianceEstimate = boost::multiprecision::sqrt(args.varianceEstimate/args.n);
		args.relativeErrorEstimate = args.sqrtVarianceEstimate / args.estimateFirstMoment;
	}
}
