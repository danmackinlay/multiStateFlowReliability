#ifndef GENERALIZED_SPLITTING_FIXED_FACTORS_HEADER_GUARD
#define GENERALIZED_SPLITTING_FIXED_FACTORS_HEADER_GUARD
#include "context.h"
namespace multistateTurnip
{
	struct generalisedSplittingFixedFactorsArgs
	{
		generalisedSplittingFixedFactorsArgs(const Context& context)
			:context(context)
		{}
		std::vector<double> levels;
		const Context& context;
		int n;
		boost::mt19937 randomSource;
		double estimate;
		std::vector<double> levelProbabilities;
		std::vector<int> splittingFactors;
		std::function<void(unsigned long, unsigned long)> progressFunction;
		std::function<void(std::string&)> outputFunc;
		double estimatedVariance;
	};
	void generalisedSplittingFixedFactors(generalisedSplittingFixedFactorsArgs& args);
}
#endif
