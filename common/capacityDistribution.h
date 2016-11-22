#ifndef CAPACITY_DISTRIBUTION_HEADER_GUARD
#define CAPACITY_DISTRIBUTION_HEADER_GUARD
#include <algorithm>
#include "includeMPFR.h"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/noncopyable.hpp>
namespace multistateTurnip
{
	class capacityDistribution : public boost::noncopyable
	{
	public:
		capacityDistribution(capacityDistribution&& other);
		capacityDistribution(const std::vector<std::pair<mpfr_class, mpfr_class>>& data);
		const std::vector<std::pair<mpfr_class, mpfr_class> >& getData() const;
		const std::vector<std::pair<double, double> >& getCumulativeData() const;
		double operator()(boost::mt19937& randomSource) const;
		double sampleConditionalLessThan(boost::mt19937& randomSource, double smaller) const;
		capacityDistribution& operator=(capacityDistribution&& other);
		capacityDistribution();
		capacityDistribution truncateAtMax(mpfr_class threshold);
		capacityDistribution createCopy();
		capacityDistribution makeCopy() const;
	private:
		capacityDistribution(const capacityDistribution& other);
		capacityDistribution& operator=(const capacityDistribution& other);
		std::vector<std::pair<mpfr_class, mpfr_class> > data;
		std::vector<std::pair<double, double> > cumulativeData;
		std::vector<std::pair<double, std::vector<std::pair<double, double> > > > conditionalCumulativeData;
	};
}
#endif
