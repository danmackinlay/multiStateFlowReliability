#include "turnip.h"
#include "computeConditionalProb.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <algorithm>
#include <boost/random/exponential_distribution.hpp>
#include "allPointsMaxFlow.hpp"
#include "edmondsKarp.h"
namespace multistateTurnip
{
	namespace turnipPrivate
	{
		struct edgeRepairData
		{
			int parallelEdgeIndex;
			double time;
			mpfr_class rate;
		};
		bool timeSorter(const edgeRepairData& first, const edgeRepairData& second)
		{
			return first.time < second.time;
		}
	}
	void turnip(turnipArgs& args)
	{
		const Context& context = args.context;
		const capacityDistribution& distribution = context.getDistribution();
		const Context::internalDirectedGraph& directedGraph = context.getDirectedGraph();
		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nVertices = boost::num_vertices(directedGraph), nDirectedEdges = boost::num_edges(directedGraph), nUndirectedEdges = boost::num_edges(graph);

		const std::vector<std::pair<double, double> >& cumulativeData = distribution.getCumulativeData();
		std::size_t nLevels = distribution.getData().size();
		std::size_t nParallelEdges = nUndirectedEdges * (cumulativeData.size()-1);
		int source = context.getSource();
		int sink = context.getSink();

		//These are only needed for the allPointsMaxFlow call. 
		allPointsMaxFlow::allPointsMaxFlowScratch<Context::internalDirectedGraph> scratch;
		//Working data for the edmonds karp call(s)
		edmondsKarpMaxFlowScratch edmondsKarpScratch;
		//All-points max-flow matrix
		std::vector<double> flowMatrix(nVertices*nVertices);
		//Get out the vertices for every undirected edge, in a vector. The index in the vector is the edge index of the edge
		std::vector<std::pair<Context::internalGraph::vertex_descriptor, Context::internalGraph::vertex_descriptor> > verticesPerEdge(nDirectedEdges);
		{
			Context::internalGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(graph);
			for(; current != end; current++)
			{
				int edgeIndex = boost::get(boost::edge_index, graph, *current);
				verticesPerEdge[edgeIndex] = std::make_pair(current->m_source, current->m_target);
			}
		}

		//Describe all the parrallel edges - In terms of the rate of that parallel edge, 
		//which original (undirected) edge it belongs to, and which index of parallel edge it is among all the parallel edges for that original edge
		std::vector<int> originalEdgeIndex;
		std::vector<int> originalEdgeLevel;
		std::vector<double> originalRates;
		std::vector<mpfr_class> originalRatesExact;
		//The initial rate at the start of each PMC step
		mpfr_class sumAllRates = 0;
		//Set up the original edge indices and rates
		for(std::size_t i = 0; i < nUndirectedEdges; i++)
		{
			mpfr_class cumulativeRates = 0;
			originalEdgeIndex.insert(originalEdgeIndex.end(), cumulativeData.size()-1, (int)i);
			originalEdgeLevel.insert(originalEdgeLevel.end(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)cumulativeData.size()-1));
			for(std::size_t j = 0; j < cumulativeData.size()-1; j++)
			{
				mpfr_class newRate = -boost::multiprecision::log(mpfr_class(cumulativeData[cumulativeData.size() - j - 1].second)) - cumulativeRates;
				originalRatesExact.push_back(newRate);
				originalRates.push_back((double)newRate);
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
		std::vector<turnipPrivate::edgeRepairData> repairTimes;
		repairTimes.reserve(nParallelEdges);
		//The edges which have already been seen. This is used to exclude edges which become redundant. 
		std::vector<bool> alreadySeen(nParallelEdges);
		//Only warn about stability once
		bool warnedStability = false;
		//The rates for all the different edges, after some edges have been discarded and their rates added to some other edge
		std::vector<mpfr_class> ratesForEdges(nParallelEdges);

		//We store the per edge repair times seperately to start with, for the purposes of stripping out stuff that's not going to be important
		std::vector<double> perEdgeRepairTimes(nLevels - 1);

		std::vector<mpfr_class> computeConditionalProbScratch;
		for(int i = 0; i < args.n; i++)
		{
			repairTimes.clear();
			std::fill(ratesForEdges.begin(), ratesForEdges.end(), 0);
			//Simulate permutation via the repair times
			for(int k = 0; k < (int)nUndirectedEdges; k++)
			{
				for(int j = 0; j < (int)nLevels - 1; j++)
				{
					boost::exponential_distribution<> repairDist(originalRates[k*(nLevels - 1) + j]);
					perEdgeRepairTimes[j] = repairDist(args.randomSource);
				}
				//The increase to highest capacity definitely occurs at some point
				turnipPrivate::edgeRepairData highest;
				highest.time = perEdgeRepairTimes[0];
				highest.rate = originalRatesExact[k*(nLevels - 1) + 0];
				highest.parallelEdgeIndex = k*(nLevels - 1) + 0;
				repairTimes.push_back(highest);

				//We can store pointers because there is no reallocation of this vector, due to reserve call
				turnipPrivate::edgeRepairData* minRepairTime = &*repairTimes.rbegin();
				for(int j = 1; j < (int)nLevels - 1; j++)
				{
					if(perEdgeRepairTimes[j] > minRepairTime->time)
					{
						minRepairTime->rate += originalRatesExact[k*(nLevels - 1) + j];
					}
					else
					{
						ratesForEdges[minRepairTime->parallelEdgeIndex] = minRepairTime->rate;
						turnipPrivate::edgeRepairData time;
						time.time = perEdgeRepairTimes[j];
						time.parallelEdgeIndex = k*(nLevels - 1) + j;
						time.rate = originalRatesExact[k*(nLevels - 1) + j];
						repairTimes.push_back(time);
						minRepairTime = &*repairTimes.rbegin();
					}
				}
				ratesForEdges[minRepairTime->parallelEdgeIndex] = minRepairTime->rate;
			}
			std::sort(repairTimes.begin(), repairTimes.end(), turnipPrivate::timeSorter);
			//No edges have yet been seen
			std::fill(alreadySeen.begin(), alreadySeen.end(), false);
			//The first rate is going to be sumAllRates
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			std::vector<turnipPrivate::edgeRepairData>::iterator repairTimeIterator = repairTimes.begin();
			//have we reached the point where we've got sufficient flow?
			bool insufficientFlow = true;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//The capacities, residual and flow are initially zero
			std::fill(capacityVector.begin(), capacityVector.end(), 0);
			std::fill(residualVector.begin(), residualVector.end(), 0);
			std::fill(flowVector.begin(), flowVector.end(), 0);
			double currentFlow = 0;
			//Counter used to make sure we only call the all-points max flow once for every fixed number of steps
			int allPointsMaxFlowCounter = 1;
			while(insufficientFlow && repairTimeIterator != repairTimes.end())
			{
				//get out the parallel edge index
				int parallelEdgeIndex = repairTimeIterator->parallelEdgeIndex;
				//Which original edge does this correspond to?
				int originalEdgeIndexThisLoop = originalEdgeIndex[parallelEdgeIndex];
				if(!alreadySeen[parallelEdgeIndex])
				{
					//Increase the capacity. Because we always remove edges with lower capacity, we don't need a maximum here
					double newCapacity = (cumulativeData.rbegin() + originalEdgeLevel[parallelEdgeIndex])->first;
					double increase = newCapacity - capacityVector[2 * originalEdgeIndexThisLoop];
					capacityVector[2 * originalEdgeIndexThisLoop] = capacityVector[2 * originalEdgeIndexThisLoop + 1] = newCapacity;
					residualVector[2 * originalEdgeIndexThisLoop] += increase;
					residualVector[2 * originalEdgeIndexThisLoop + 1] += increase;

					//determine whether or not we've hit the critical threshold
					edmondsKarpMaxFlow(&capacityVector.front(), &flowVector.front(), &residualVector.front(), directedGraph, source, sink, args.threshold, edmondsKarpScratch, currentFlow);
					insufficientFlow = args.threshold > currentFlow;
					//Add the current rate
					ratesForPMC.push_back(currentRate);
					currentRate -= ratesForEdges[parallelEdgeIndex];
					alreadySeen[parallelEdgeIndex] = true;
					//If we now have newThreshold or higher flow between the the vertices for edge originalEdgeIndexThisLoop,
					//we can discard ALL not yet added edges between those two vertices
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
											currentRate -= ratesForEdges[parallelEdgeIndex];
											alreadySeen[parallelEdgeIndex] = true;
										}
										parallelEdgeIndex++;
									} while (parallelEdgeIndex >= 0 && parallelEdgeIndex < (int)nParallelEdges);
								}
							}
						}
					}
					else 
					{
						Context::internalGraph::vertex_descriptor firstVertex = verticesPerEdge[originalEdgeIndexThisLoop].first;
						Context::internalGraph::vertex_descriptor secondVertex = verticesPerEdge[originalEdgeIndexThisLoop].second;
						//We want to start this computation from the beginning (rather than the incremental version). So the flows start as zero
						std::fill(flowVectorIncreasedEdge.begin(), flowVectorIncreasedEdge.end(), 0);
						std::copy(capacityVector.begin(), capacityVector.end(), residualVectorIncreasedEdge.begin());
						double flowIncreasedEdge = 0;
						edmondsKarpMaxFlow(&capacityVector.front(), &flowVectorIncreasedEdge.front(), &residualVectorIncreasedEdge.front(), directedGraph, firstVertex, secondVertex, args.threshold, edmondsKarpScratch, flowIncreasedEdge);
						if(flowIncreasedEdge >= args.threshold)
						{
							parallelEdgeIndex = (int)(originalEdgeIndexThisLoop*(nLevels-1));
							do
							{
								if(originalEdgeIndex[parallelEdgeIndex] != originalEdgeIndexThisLoop) break;
								if(!alreadySeen[parallelEdgeIndex]) 
								{
									currentRate -= ratesForEdges[parallelEdgeIndex];
									alreadySeen[parallelEdgeIndex] = true;
								}
								parallelEdgeIndex++;
							}
							while(parallelEdgeIndex >= 0 && parallelEdgeIndex < (int)nParallelEdges);
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
