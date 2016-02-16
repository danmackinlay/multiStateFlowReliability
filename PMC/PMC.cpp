#include <boost/program_options.hpp>
#include "computeConditionalProb.h"
#include "arguments.h"
#include "argumentsMPFR.h"
#include "includeMPFR.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <algorithm>
#include <boost/random/exponential_distribution.hpp>
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
			("turnip", "(Flag) Should we make turnip-style optimisations?")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
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
		bool useTurnip = variableMap.count("turnip") > 0;
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
		Context context = Context::emptyContext();
		double newThresholdD;
		{
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
			newThresholdD = (double)newThreshold;
			//truncate the set of possible capacities at the new threshold, otherwise
			//sorting that many edge repair times is a bottleneck in itself. 
			movedDistribution = movedDistribution.truncateAtMax(newThreshold);
			if(newThreshold < 0)
			{
				std::cout << "Minimum flow for all links is higher than the threshold" << std::endl;
				return 0;
			}
			if(!readContext(variableMap, context, std::move(movedDistribution), newThreshold))
			{
				return 0;
			}
		}

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);
		
		const capacityDistribution& distribution = context.getDistribution();
		const std::vector<std::pair<double, double> >& cumulativeData = distribution.getCumulativeData();
		std::size_t nLevels = distribution.getData().size();
		std::size_t nParallelEdges = context.getNEdges() * (nLevels-1);

		//Describe all the parrallel edges - In terms of the rate of that parallel edge, 
		//which original edge it belongs to, and which index of parallel edge it is among all the parallel edges 
		//for that original edge
		std::vector<int> originalEdgeIndex;
		std::vector<int> originalEdgeLevel;
		std::vector<double> originalRates;
		std::vector<mpfr_class> originalRatesExact;
		//Set up the original edge indices and rates
		//The initial rate at the start of each PMC step
		mpfr_class sumAllRates = 0;
		for(std::size_t i = 0; i < context.getNEdges(); i++)
		{
			mpfr_class cumulativeRates = 0;
			originalEdgeIndex.insert(originalEdgeIndex.end(), nLevels-1, (int)i);
			originalEdgeLevel.insert(originalEdgeLevel.end(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nLevels-1));
			for(std::size_t j = 0; j < nLevels-1; j++)
			{
				mpfr_class newRate = -boost::multiprecision::log(mpfr_class(cumulativeData[nLevels - j - 1].second)) - cumulativeRates;
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
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC;
		//Repair time vector
		std::vector<std::pair<int, double> > repairTimes;
		repairTimes.resize(nParallelEdges);
		//The edges which have already been seen. This is used to exclude edges which become redundant. 
		std::vector<bool> alreadySeen(nParallelEdges);
		//Only warn about stability once
		bool warnedStability = false;
		std::vector<mpfr_class> computeConditionalProbScratch;
		for(int i = 0; i < n; i++)
		{
			//Simulate permutation
			for(int j = 0; j < (int)nParallelEdges; j++)
			{
				boost::exponential_distribution<> repairDist(originalRates[j]);
				repairTimes[j].second = repairDist(randomSource);
				repairTimes[j].first = j;
			}
			std::sort(repairTimes.begin(), repairTimes.end(), secondArgumentSorter);
			//No edges have yet been seen
			std::fill(alreadySeen.begin(), alreadySeen.end(), false);
			//The first rate is going to be this
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			int permutationCounter = 0;
			//have we reached the point where we've got sufficient flow?
			bool insufficientFlow = true;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//The capacities are initially zero
			std::fill(capacityVector.begin(), capacityVector.end(), 0);
			while(insufficientFlow)
			{
				//get out the parallel edge index
				int parallelEdgeIndex = repairTimes[permutationCounter].first;
				//Which original edge does this correspond to?
				int originalEdgeIndexThisLoop = originalEdgeIndex[parallelEdgeIndex];
				//If we're using the partial turnip option this edge might have been discounted already
				//because it doesn't add anything to the maximum flow
				if(!alreadySeen[parallelEdgeIndex])
				{
					//Increase the capacity
					capacityVector[2 * originalEdgeIndexThisLoop] = capacityVector[2 * originalEdgeIndexThisLoop + 1] = std::max(capacityVector[2 * originalEdgeIndexThisLoop + 1], (cumulativeData.rbegin() + originalEdgeLevel[parallelEdgeIndex])->first);
					//determine whether or not we've hit the critical threshold
					double currentFlow = context.getMaxFlow(capacityVector);
					insufficientFlow = newThresholdD > currentFlow;
					//Add the current rate
					ratesForPMC.push_back(currentRate);
					if(useTurnip)
					{
						//Start going back through the other parallel edges for this original edge
						//if they hadn't already been observed to occur, adjust the rate to indicate that we don't need them.
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
					else
					{
						currentRate -= originalRates[parallelEdgeIndex];
						alreadySeen[parallelEdgeIndex] = true;
					}
				}
				permutationCounter++;
			}
			mpfr_class additionalPart = computeConditionalProb(ratesForPMC, computeConditionalProbScratch);
			//mpfr_class additionalPart2 = computeConditionalProb(ratesForPMC);
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
