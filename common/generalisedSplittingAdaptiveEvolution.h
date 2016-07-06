#ifndef GENERALIZED_SPLITTING_ADAPTIVE_EVOLUTION_HEADER_GUARD
#define GENERALIZED_SPLITTING_ADAPTIVE_EVOLUTION_HEADER_GUARD
#include "context.h"
namespace multistateTurnip
{
	struct generalisedSplittingAdaptiveEvolutionArgs
	{
		generalisedSplittingAdaptiveEvolutionArgs(const Context& context)
			:context(context)
		{}
		std::vector<double> times;
		double level;
		const Context& context;
		int n;
		int fraction;
		boost::mt19937 randomSource;
		std::function<void(unsigned long, unsigned long)> progressFunction;
		std::function<void(std::string&)> outputFunc;
	};
	void generalisedSplittingAdaptiveEvolution(generalisedSplittingAdaptiveEvolutionArgs& args);
}
#endif
