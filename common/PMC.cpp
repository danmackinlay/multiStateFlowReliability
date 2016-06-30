#include "PMC.h"
#include "computeConditionalProb.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <algorithm>
#include <boost/random/exponential_distribution.hpp>
#include "edmondsKarp.hpp"
#include "turnip.h"
namespace multistateTurnip
{
	void pmc(pmcArgs& args)
	{
		const Context::internalDirectedGraph& directedGraph = args.context.getDirectedGraph();
		const Context::internalGraph& undirectedGraph = args.context.getGraph();

		std::size_t nDirectedEdges = boost::num_edges(directedGraph), nUndirectedEdges = boost::num_edges(undirectedGraph);
		int source = args.context.getSource();
		int sink = args.context.getSink();

		//Working memory for edmonds Karp
		edmondsKarpMaxFlowScratch<Context::internalDirectedGraph, double> scratch;

		std::vector<double> minimumCapacities(nUndirectedEdges*2);

		//Rates for the parallel edges. Different parallel edges for the same underlying original edge are consecutive
		std::vector<std::vector<double> > originalRates(nUndirectedEdges);
		std::vector<std::vector<mpfr_class> > originalRatesExact(nUndirectedEdges);
		//The initial rate at the start of each PMC step
		mpfr_class sumAllRates = 0;
		//Set up the edge rates
		int totalLevels = 0;
		for(std::size_t i = 0; i < nUndirectedEdges; i++)
		{
			std::vector<double>& currentEdgeRates = originalRates[i];
			std::vector<mpfr_class>& currentEdgeRatesExact = originalRatesExact[i];
			const capacityDistribution& currentEdgeDistribution = args.context.getDistribution(i);
			const std::vector<std::pair<double, double> >& cumulativeData = currentEdgeDistribution.getCumulativeData();
			std::size_t nLevels = currentEdgeDistribution.getData().size();
			minimumCapacities[2*i] = minimumCapacities[2*i+1] = cumulativeData.front().first;
			totalLevels += nLevels-1;
			mpfr_class cumulativeRates = 0;
			for(std::size_t j = 0; j < nLevels-1; j++)
			{
				mpfr_class newRate = -boost::multiprecision::log(mpfr_class(cumulativeData[nLevels - j - 1].second)) - cumulativeRates;
				currentEdgeRatesExact.push_back(newRate);
				currentEdgeRates.push_back((double)newRate);
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
		std::vector<edgeRepairData> repairTimes;
		repairTimes.reserve(totalLevels);

		//We store the per edge repair times seperately to start with, for the purposes of stripping out stuff that's not going to be important
		std::vector<double> perEdgeRepairTimes;

		//Only warn about stability once
		bool warnedStability = false;
		std::vector<mpfr_class> computeConditionalProbScratch;
		//Work out the minimum possible flow
		double minimumPossibleFlow = 0;
		std::copy(minimumCapacities.begin(), minimumCapacities.end(), capacityVector.begin());
		std::copy(minimumCapacities.begin(), minimumCapacities.end(), residualVector.begin());
		std::fill(flowVector.begin(), flowVector.end(), 0);
		edmondsKarpMaxFlow<Context::internalDirectedGraph, double>(&capacityVector.front(), &flowVector.front(), &residualVector.front(), directedGraph, source, sink, args.threshold, scratch, minimumPossibleFlow);
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
			for(int k = 0; k < (int)nUndirectedEdges; k++)
			{
				std::vector<double>& currentEdgeRates = originalRates[k];
				std::vector<mpfr_class>& currentEdgeRatesExact = originalRatesExact[k];
				const capacityDistribution& currentEdgeDistribution = args.context.getDistribution(k);
				std::size_t nLevels = currentEdgeDistribution.getData().size();
				perEdgeRepairTimes.resize(nLevels-1);
				//Simulate the times for all the parallel edges, for a single underlying edge
				for(int j = 0; j < (int)nLevels - 1; j++)
				{
					boost::exponential_distribution<> repairDist(currentEdgeRates[j]);
					perEdgeRepairTimes[j] = repairDist(args.randomSource);
				}
				//The increase to highest capacity definitely occurs at some point
				edgeRepairData highest;
				highest.time = perEdgeRepairTimes[0];
				highest.rate = currentEdgeRates[0];
				highest.edge = k;
				highest.level = 0;
				repairTimes.push_back(highest);

				//We can store pointers because there is no reallocation of this vector, due to reserve call
				edgeRepairData* minRepairTime = &repairTimes.back();
				//Aggregate rates. If a lower capacity edge occurs later, add its rate to the current minimum repair time edge. 
				for(int j = 1; j < (int)nLevels - 1; j++)
				{
					if(perEdgeRepairTimes[j] > minRepairTime->time)
					{
						minRepairTime->rate += currentEdgeRatesExact[j];
					}
					else
					{
						edgeRepairData time;
						time.time = perEdgeRepairTimes[j];
						time.level = j;
						time.rate = currentEdgeRatesExact[j];
						time.edge = k;
						repairTimes.push_back(time);
						minRepairTime = &repairTimes.back();
					}
				}
			}
			std::sort(repairTimes.begin(), repairTimes.end(), timeSorter);
			//The first rate is going to be this
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
			double currentFlow = minimumPossibleFlow;
			while(insufficientFlow && repairTimeIterator != repairTimes.end())
			{
				//get out the parallel edge index
				int level = repairTimeIterator->level;
				//Which original edge does this correspond to?
				int edge = repairTimeIterator->edge;
				const capacityDistribution& relevantDistribution = args.context.getDistribution(edge);
				const std::vector<std::pair<double, double> >& cumulativeData = relevantDistribution.getCumulativeData();
				//Increase the capacity
				double newCapacity = std::max(capacityVector[2 * edge + 1], (cumulativeData.rbegin() + level)->first);
				double increase = newCapacity - capacityVector[2 * edge];
				capacityVector[2 * edge] = capacityVector[2 * edge + 1] = newCapacity;
				residualVector[2 * edge] += increase;
				residualVector[2 * edge + 1] += increase;
				//determine whether or not we've hit the critical threshold
				edmondsKarpMaxFlow<Context::internalDirectedGraph, double>(&capacityVector.front(), &flowVector.front(), &residualVector.front(), directedGraph, source, sink, args.threshold, scratch, currentFlow);
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
