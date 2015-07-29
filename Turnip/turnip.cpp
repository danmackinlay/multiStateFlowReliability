#include <boost/program_options.hpp>
#include "computeConditionalProb.h"
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "includeMPFR.h"
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
	bool secondArgumentSorter(const std::pair<int, double>& first, const std::pair<int, double>& second)
	{
		return first.second < second.second;
	}
	int main(int argc, char **argv)
	{
		mpfr_set_default_prec(1024);

		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and completeGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph and completeGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("capacityFile", boost::program_options::value<std::string>(), "(string) File containing capacity distribution data. Incompatible with graphFile and gridGraph")
			("uniformCapacity", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) Pair of integers a and b, so that capacities are uniform on the integer range [a, b]")
			("threshold", boost::program_options::value<std::string>(), "(float) The threshold flow between the source and sink")
			("n", boost::program_options::value<int>(), "(int) The number of simulations to perform. ")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("full", boost::program_options::value<int>()->implicit_value(1), "(int) Run all-points max flow once for every n added edges. Defaults to 1. ")
			("help", "Display this message");
		
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << 
				"This program estimates the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component?\n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}

		int n;
		if(!readN(variableMap, n))
		{
			return 0;
		}
		mpfr_class threshold;
		if(!readThreshold(variableMap, threshold))
		{
			std::cout << "Unable to parse threshold value" << std::endl;
			return 0;
		}
		//Don't use this distribution - It gets moved into the Context object
		capacityDistribution movedDistribution;
		std::string error;
		mpfr_class originalMinFlow;
		//This alter the minimum flow to be zero, and returns the original.
		if(!readCapacityDistribution(variableMap, movedDistribution, originalMinFlow, error))
		{
			std::cout << error << std::endl;
			return 0;
		}
		mpfr_class newThreshold = threshold - originalMinFlow;
		double newThresholdD = (double)newThreshold;
		//truncate the set of possible capacities at the new threshold, otherwise
		//sorting that many edge repair times is a bottleneck in itself. 
		movedDistribution = movedDistribution.truncateAtMax(newThreshold);
		if(newThreshold < 0)
		{
			std::cout << "Minimum flow for all links is higher than the threshold" << std::endl;
			return 0;
		}
		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, movedDistribution, newThreshold))
		{
			return 0;
		}
		bool fullSpecified = variableMap.count("full") > 0;
		int fullAllPointsIncrement = variableMap["full"].as<int>();

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);
		
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
		for(int i = 0; i < n; i++)
		{
			//Simulate permutation
			for(int j = 0; j < nParallelEdges; j++)
			{
				boost::exponential_distribution<> repairDist(originalRates[j]);
				repairTimes[j].second = repairDist(randomSource);
				repairTimes[j].first = j;
			}
			std::sort(repairTimes.begin(), repairTimes.end(), secondArgumentSorter);
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
					insufficientFlow = newThresholdD > currentFlow;
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
					} while (parallelEdgeIndex >= 0 && parallelEdgeIndex < nParallelEdges);

					//If we now have newThreshold or higher flow between the the vertices for edge originalEdgeIndexThisLoop,
					//we can discard ALL not yet added edges between those two vertices
					Context::internalGraph::vertex_descriptor firstVertex = verticesPerEdge[originalEdgeIndexThisLoop].first;
					Context::internalGraph::vertex_descriptor secondVertex = verticesPerEdge[originalEdgeIndexThisLoop].second;
					if (fullSpecified)
					{
						if ((allPointsMaxFlowCounter++ % fullAllPointsIncrement) == 0)
						{
							allPointsMaxFlow::allPointsMaxFlow<Context::internalDirectedGraph>(flowMatrix, capacityVector, directedGraph, scratch);
							Context::internalGraph::edge_iterator current, end;
							boost::tie(current, end) = boost::edges(graph);
							for (; current != end; current++)
							{
								Context::internalGraph::vertex_descriptor source = current->m_source, target = current->m_target;
								if (flowMatrix[source + nVertices * target] >= newThresholdD)
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
									} while (parallelEdgeIndex >= 0 && parallelEdgeIndex < nParallelEdges);
								}
							}
						}
					}
					else if(context.getMaxFlow(capacityVector, firstVertex, secondVertex) >= newThresholdD)
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
						while(parallelEdgeIndex >= 0 && parallelEdgeIndex < nParallelEdges);
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
		mpfr_class estimateFirstMoment = sumConditional/n;
		mpfr_class estimateSecondMoment = sumSquaredConditional/n;
		mpfr_class varianceEstimate = estimateSecondMoment - estimateFirstMoment*estimateFirstMoment;
		mpfr_class sqrtVarianceEstimate = boost::multiprecision::sqrt(varianceEstimate/n);
		mpfr_class relativeErrorEstimate = sqrtVarianceEstimate / estimateFirstMoment;

		std::cout << "Unreliability probability estimate was " << (double)estimateFirstMoment << std::endl;
		std::cout << "Relative error was " << (double)relativeErrorEstimate << std::endl;
		return 0;
	}
}

int main(int argc, char **argv)
{
	return multistateTurnip::main(argc, argv);
}