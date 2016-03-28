#ifndef GENERALIZED_SPLITTING_HEADER_GUARD
#define GENERALIZED_SPLITTING_HEADER_GUARD
#include "context.h"
namespace multistateTurnip
{
	struct generalisedSplittingArgs
	{
		generalisedSplittingArgs(const Context& context)
			:context(context)
		{}
		std::vector<double> levels;
		const Context& context;
		int n;
		boost::mt19937 randomSource;
		double estimate;
		std::function<void(unsigned long, unsigned long)> progressFunction;
		std::function<void(std::string&)> outputFunc;
	};
	void generalisedSplitting(generalisedSplittingArgs& args);
}
#endif