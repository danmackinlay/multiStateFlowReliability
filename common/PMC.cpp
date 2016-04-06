#include "PMC.h"
#include "computeConditionalProb.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <algorithm>
#include <boost/random/exponential_distribution.hpp>
#include "edmondsKarp.h"
namespace multistateTurnip
{
	namespace pmcPrivate
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
	void pmc(pmcArgs& args)
	{
		const capacityDistribution& distribution = args.context.getDistribution();
		const std::vector<std::pair<double, double> >& cumulativeData = distribution.getCumulativeData();
		const Context::internalDirectedGraph& directedGraph = args.context.getDirectedGraph();
		const Context::internalGraph& undirectedGraph = args.context.getGraph();

		std::size_t nDirectedEdges = boost::num_edges(directedGraph), nUndirectedEdges = boost::num_edges(undirectedGraph);
		std::size_t nLevels = distribution.getData().size();
		std::size_t nParallelEdges = nUndirectedEdges * (nLevels-1);
		int source = args.context.getSource();
		int sink = args.context.getSink();

		//Working memory for edmonds Karp
		edmondsKarpMaxFlowScratch scratch;

		//Describe all the parallel edges - In terms of the rate of that parallel edge, 
		//which original (undirected) edge it belongs to, and which index of parallel edge it is among all the parallel edges for that original edge
		std::vector<int> originalEdgeIndex;
		std::vector<int> originalEdgeLevel;
		//Rates for the parallel edges. Different parallel edges for the same underlying original edge are consecutive
		std::vector<double> originalRates;
		std::vector<mpfr_class> originalRatesExact;
		//Set up the original edge indices and rates
		//The initial rate at the start of each PMC step
		mpfr_class sumAllRates = 0;
		for(std::size_t i = 0; i < nUndirectedEdges; i++)
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
		//Residual and flow vectors for edmonds-karp
		std::vector<double> residualVector(nDirectedEdges);
		std::vector<double> flowVector(nDirectedEdges);
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC;
		//Repair time vector
		std::vector<pmcPrivate::edgeRepairData> repairTimes;
		repairTimes.reserve(nParallelEdges);

		//We store the per edge repair times seperately to start with, for the purposes of stripping out stuff that's not going to be important
		std::vector<double> perEdgeRepairTimes(nLevels - 1);

		//Only warn about stability once
		bool warnedStability = false;
		std::vector<mpfr_class> computeConditionalProbScratch;
		for(int i = 0; i < args.n; i++)
		{
			repairTimes.clear();
			double currentFlow = 0;
			//Simulate permutation via the repair times
			for(int k = 0; k < (int)nUndirectedEdges; k++)
			{
				//Simulate the times for all the parallel edges, for a single underlying edge
				for(int j = 0; j < (int)nLevels - 1; j++)
				{
					boost::exponential_distribution<> repairDist(originalRates[k*(nLevels - 1) + j]);
					perEdgeRepairTimes[j] = repairDist(args.randomSource);
				}
				//The increase to highest capacity definitely occurs at some point
				pmcPrivate::edgeRepairData highest;
				highest.time = perEdgeRepairTimes[0];
				highest.rate = originalRatesExact[k*(nLevels - 1) + 0];
				highest.parallelEdgeIndex = k*(nLevels - 1) + 0;
				repairTimes.push_back(highest);

				//We can store pointers because there is no reallocation of this vector, due to reserve call
				pmcPrivate::edgeRepairData* minRepairTime = &repairTimes.back();
				//Aggregate rates. If a lower capacity edge occurs later, add its rate to the current minimum repair time edge. 
				for(int j = 1; j < (int)nLevels - 1; j++)
				{
					if(perEdgeRepairTimes[j] > minRepairTime->time)
					{
						minRepairTime->rate += originalRatesExact[k*(nLevels - 1) + j];
					}
					else
					{
						pmcPrivate::edgeRepairData time;
						time.time = perEdgeRepairTimes[j];
						time.parallelEdgeIndex = k*(nLevels - 1) + j;
						time.rate = originalRatesExact[k*(nLevels - 1) + j];
						repairTimes.push_back(time);
						minRepairTime = &repairTimes.back();
					}
				}
			}
			std::sort(repairTimes.begin(), repairTimes.end(), pmcPrivate::timeSorter);
			//The first rate is going to be this
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			std::vector<pmcPrivate::edgeRepairData>::iterator repairTimeIterator = repairTimes.begin();
			//have we reached the point where we've got sufficient flow?
			bool insufficientFlow = true;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//The capacities, residuals and flows are initially zero
			std::fill(capacityVector.begin(), capacityVector.end(), 0);
			std::fill(residualVector.begin(), residualVector.end(), 0);
			std::fill(flowVector.begin(), flowVector.end(), 0);
			while(insufficientFlow && repairTimeIterator != repairTimes.end())
			{
				//get out the parallel edge index
				int parallelEdgeIndex = repairTimeIterator->parallelEdgeIndex;
				//Which original edge does this correspond to?
				int originalEdgeIndexThisLoop = originalEdgeIndex[parallelEdgeIndex];
				//Increase the capacity
				double newCapacity = std::max(capacityVector[2 * originalEdgeIndexThisLoop + 1], (cumulativeData.rbegin() + originalEdgeLevel[parallelEdgeIndex])->first);
				double increase = newCapacity - capacityVector[2 * originalEdgeIndexThisLoop];
				capacityVector[2 * originalEdgeIndexThisLoop] = capacityVector[2 * originalEdgeIndexThisLoop + 1] = newCapacity;
				residualVector[2 * originalEdgeIndexThisLoop] += increase;
				residualVector[2 * originalEdgeIndexThisLoop + 1] += increase;
				//determine whether or not we've hit the critical threshold
				edmondsKarpMaxFlow(&capacityVector.front(), &flowVector.front(), &residualVector.front(), directedGraph, source, sink, args.threshold, scratch, currentFlow);
				insufficientFlow = args.threshold > currentFlow;
				//Add the current rate
				ratesForPMC.push_back(currentRate);
				currentRate -= repairTimeIterator->rate;
				repairTimeIterator++;
			}
			mpfr_class additionalPart;
			//compute conditional probability
			if(!insufficientFlow)
			{
				additionalPart = computeConditionalProb(ratesForPMC, computeConditionalProbScratch);
				//mpfr_class additionalPart2 = computeConditionalProb(ratesForPMC);
				if(additionalPart > 1 && !warnedStability)
				{
					std::string output = "Numerical stability problem detected";
					args.outputFunc(output);
					warnedStability = true;
				}
			}
			else additionalPart = 1;
			//Add this conditional probability to the running sum and running sum of squares
			sumConditional += additionalPart;
			sumSquaredConditional += additionalPart*additionalPart;
			if(i % 100 == 0)
			{
				args.progressFunction(i, args.n);
			}
		}
		//Return average of all the terms in the sum
		args.firstMomentSingleSample = sumConditional/args.n;
		//Return the estimated second moment of the terms in the sum
		args.secondMomentSingleSample = sumSquaredConditional/args.n;
		//The variance of a single term in the sum
		args.varianceSingleSample = args.secondMomentSingleSample - args.firstMomentSingleSample*args.firstMomentSingleSample;
		args.sqrtVarianceOfEstimate = boost::multiprecision::sqrt(args.varianceSingleSample/args.n);
		args.relativeErrorEstimate = args.sqrtVarianceOfEstimate / args.firstMomentSingleSample;
	}
}
