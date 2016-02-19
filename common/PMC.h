#ifndef PMC_HEADER_GUARD
#define PMC_HEADER_GUARD
#include <functional>
#include "context.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
namespace multistateTurnip
{
	struct pmcArgs
	{
		pmcArgs(const Context& context)
			:context(context)
		{}
		const Context& context;
		mpfr_class estimateFirstMoment, estimateSecondMoment, varianceEstimate, sqrtVarianceEstimate, relativeErrorEstimate;
		int n;
		double threshold;
		std::function<void(std::string&)> outputFunc;
		boost::mt19937 randomSource;
	};
	void pmc(pmcArgs& args);
}
#endif
