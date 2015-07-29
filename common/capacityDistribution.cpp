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
		cumulativeData.resize(data.size());
		mpfr_class sum = 0;
		for(std::size_t i = 0; i < data.size(); i++)
		{
			cumulativeData[i].first = (double)data[i].first;
			cumulativeData[i].second = (double)sum;
			sum += data[i].second;
		}
		if(cumulativeData.size() > 1 && (cumulativeData.rbegin()+1)->second > 1)
		{
			throw std::runtime_error("The sum of the first n - 1 probabilities was already bigger than 1");
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
		return *this;
	}
	capacityDistribution::capacityDistribution(capacityDistribution&& other)
	{
		data.swap(other.data);
		cumulativeData.swap(other.cumulativeData);
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
		truncatedData.push_back(std::make_pair(newThreshold, thresholdProbability));
		return capacityDistribution(truncatedData);
	}
}