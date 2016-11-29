#ifndef TURNIP_DIRECTED_HEADER_GUARD
#define TURNIP_DIRECTED_HEADER_GUARD
#include "contextDirected.h"
#include <boost/random/mersenne_twister.hpp>
#include <functional>
#include "includeMPFR.h"
#include "turnip.h"
namespace multistateTurnip
{
	struct turnipDirectedArgs
	{
	public:
		turnipDirectedArgs(const ContextDirected& context)
			:context(context)
		{}
		const ContextDirected& context;
		int n;
		double threshold;
		std::function<void(std::string&)> outputFunc;
		boost::mt19937 randomSource;
		mpfr_class firstMomentSingleSample, secondMomentSingleSample, varianceSingleSample, sqrtVarianceOfEstimate, relativeErrorEstimate;
	};
	void turnipDirected(turnipDirectedArgs& args);
}
#endif
