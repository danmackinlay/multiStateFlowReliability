#include "PMC.h"
#include "computeConditionalProb.h"
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
	void pmc(pmcArgs& args)
	{
		const capacityDistribution& distribution = args.context.getDistribution();
		const std::vector<std::pair<double, double> >& cumulativeData = distribution.getCumulativeData();
		std::size_t nLevels = distribution.getData().size();
		std::size_t nParallelEdges = args.context.getNEdges() * (nLevels-1);

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
		for(std::size_t i = 0; i < args.context.getNEdges(); i++)
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
		std::vector<double>& capacityVector = args.context.getCapacityVector();
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
		for(int i = 0; i < args.n; i++)
		{
			//Simulate permutation
			for(int j = 0; j < (int)nParallelEdges; j++)
			{
				boost::exponential_distribution<> repairDist(originalRates[j]);
				repairTimes[j].second = repairDist(args.randomSource);
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
					double currentFlow = args.context.getMaxFlow(capacityVector);
					insufficientFlow = args.threshold > currentFlow;
					//Add the current rate
					ratesForPMC.push_back(currentRate);
					if(args.useTurnip)
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
				std::string output = "Numerical stability problem detected";
				args.outputFunc(output);
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
