#ifndef TURNIP_HEADER_GUARD
#define TURNIP_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
#include <functional>
#include "includeMPFR.h"
namespace multistateTurnip
{
	struct turnipArgs
	{
	public:
		turnipArgs(const Context& context)
			:context(context)
		{}
		const Context& context;
		int n;
		double threshold;
		std::function<void(std::string&)> outputFunc;
		boost::mt19937 randomSource;
		bool useAllPointsMaxFlow;
		int allPointsMaxFlowIncrement;
		mpfr_class estimateFirstMoment, estimateSecondMoment, varianceEstimate, sqrtVarianceEstimate, relativeErrorEstimate;
	};
	void turnip(turnipArgs& args);
}
#endif
