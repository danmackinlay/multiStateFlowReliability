#ifndef PMC_DIRECTED_HEADER_GUARD
#define PMC_DIRECTED_HEADER_GUARD
#include <functional>
#include "contextDirected.h"
#include "includeMPFR.h"
#include <boost/random/mersenne_twister.hpp>
namespace multistateTurnip
{
	struct pmcDirectedArgs
	{
		pmcDirectedArgs(const ContextDirected& context)
			:context(context)
		{}
		const ContextDirected& context;
		mpfr_class firstMomentSingleSample, secondMomentSingleSample, varianceSingleSample, sqrtVarianceOfEstimate, relativeErrorEstimate;
		int n;
		double threshold;
		std::function<void(std::string&)> outputFunc;
		boost::mt19937 randomSource;
		std::function<void(unsigned long, unsigned long)> progressFunction;
	};
	void pmcDirected(pmcDirectedArgs& args);
}
#endif
