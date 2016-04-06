#include "capacityDistribution.h"
#include <algorithm>
#include <boost/random/uniform_01.hpp>
#include <boost/tuple/tuple.hpp>
namespace multistateTurnip
{
	capacityDistribution::capacityDistribution(const std::vector<std::pair<mpfr_class, mpfr_class>>& data)
		:data(data)
	{
		if(data.size() == 0)
		{
			throw std::runtime_error("Input vector cannot be empty");
		}
		//Cumulative data
		cumulativeData.resize(data.size());
		mpfr_class sum = 0;
		for(std::size_t i = 0; i < data.size(); i++)
		{
			cumulativeData[i].first = (double)data[i].first;
			cumulativeData[i].second = (double)sum;
			sum += data[i].second;
		}
		double discrepancy = (double)(mpfr_class)(sum - 1);
		if(fabs(discrepancy) > 1e-6)
		{
			throw std::runtime_error("Capacity probabilities must sum to 1");
		}
		if(cumulativeData.size() > 1 && (cumulativeData.rbegin()+1)->second > 1)
		{
			throw std::runtime_error("The sum of the first n - 1 probabilities was already bigger than 1");
		}

		//The conditional cumulative data. The last set of values calculated here are actually the same as the values in cumulativeData
		conditionalCumulativeData.resize(data.size());
		for(std::size_t i = 1; i <= data.size(); i++)
		{
			sum = 0;
			//Entry j of this vector is the probability of being less than or equal to level j, conditional on knowing that the value is less than or equal to level i-1. 
			std::vector<std::pair<double, double> > currentCumulativeData(i);
			for(std::size_t j = 0; j < i; j++)
			{
				currentCumulativeData[j].first = (double)data[j].first;
				currentCumulativeData[j].second = (double)sum;
				sum += data[j].second;
			}
			for(std::size_t j = 0; j < i; j++)
			{
				currentCumulativeData[j].second /= (double)sum;
			}
			conditionalCumulativeData[i-1].first = (double)data[i-1].first;
			conditionalCumulativeData[i-1].second.swap(currentCumulativeData);
		}
	}
	capacityDistribution::capacityDistribution()
	{}
	const std::vector<std::pair<mpfr_class, mpfr_class> >& capacityDistribution::getData() const
	{
		return data;
	}
	const std::vector<std::pair<double, double> >& capacityDistribution::getCumulativeData() const
	{
		return cumulativeData;
	}
	bool sortSecond(const std::pair<double, double>& first, const std::pair<double, double>& second)
	{
		return first.second < second.second;
	}
	double capacityDistribution::operator()(boost::mt19937& randomSource) const
	{
		boost::uniform_01<> uniform;
		double cdfValue = uniform(randomSource);

		std::vector<std::pair<double, double> >::const_iterator lowerBound, upperBound;
		std::pair<double, double> searchFor(0, cdfValue);
		boost::tie(lowerBound, upperBound) = std::equal_range(cumulativeData.begin(), cumulativeData.end(), searchFor, sortSecond);
		double value;
		if(lowerBound == cumulativeData.end()) value = cumulativeData.rbegin()->first;
		else value = (lowerBound-1)->first;
		return value;
	}
	capacityDistribution& capacityDistribution::operator=(capacityDistribution&& other)
	{
		data.swap(other.data);
		cumulativeData.swap(other.cumulativeData);
		conditionalCumulativeData.swap(other.conditionalCumulativeData);
		return *this;
	}
	capacityDistribution::capacityDistribution(capacityDistribution&& other)
	{
		data.swap(other.data);
		cumulativeData.swap(other.cumulativeData);
		conditionalCumulativeData.swap(other.conditionalCumulativeData);
	}
	capacityDistribution capacityDistribution::truncateAtMax(mpfr_class newThreshold)
	{
		mpfr_class thresholdProbability = 0;
		std::vector<std::pair<mpfr_class, mpfr_class> > truncatedData;
		for (std::vector<std::pair<mpfr_class, mpfr_class> >::const_iterator i = data.begin(); i != data.end(); i++)
		{
			if (i->first < newThreshold)
			{
				truncatedData.push_back(*i);
			}
			else thresholdProbability += i->second;
		}
		if(thresholdProbability > 0) truncatedData.push_back(std::make_pair(newThreshold, thresholdProbability));
		return capacityDistribution(truncatedData);
	}
	double capacityDistribution::sampleConditionalLessThan(boost::mt19937& randomSource, double smaller) const
	{
		boost::uniform_01<> uniform;
		double cdfValue = uniform(randomSource);
		std::pair<double, double> searchFor(0, cdfValue);
		for(std::vector<std::pair<double, std::vector<std::pair<double, double> > > >::const_reverse_iterator i = conditionalCumulativeData.rbegin(); i != conditionalCumulativeData.rend(); i++)
		{
			if(i->first < smaller)
			{
				const std::vector<std::pair<double, double> >& relevantVector = i->second;
				std::vector<std::pair<double, double> >::const_iterator lowerBound, upperBound;
				boost::tie(lowerBound, upperBound) = std::equal_range(relevantVector.begin(), relevantVector.end(), searchFor, sortSecond);
				double value;
				if(lowerBound == cumulativeData.end()) value = relevantVector.rbegin()->first;
				else if(lowerBound == cumulativeData.begin()) value = lowerBound->first;
				else value = (lowerBound-1)->first;
				return value;
			}
		}
		throw std::runtime_error("Internal error");
	}
}
