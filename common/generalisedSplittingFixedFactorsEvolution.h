#ifndef GENERALIZED_SPLITTING_FIXED_FACTORS_EVOLUTION_HEADER_GUARD
#define GENERALIZED_SPLITTING_FIXED_FACTORS_EVOLUTION_HEADER_GUARD
#include "context.h"
namespace multistateTurnip
{
	struct generalisedSplittingFixedFactorsEvolutionArgs
	{
		generalisedSplittingFixedFactorsEvolutionArgs(const Context& context)
			:context(context)
		{}
		std::vector<double> times;
		double level;
		const Context& context;
		int n;
		boost::mt19937 randomSource;
		double estimate;
		std::vector<double> timeProbabilities;
		std::vector<int> splittingFactors;
		std::function<void(unsigned long, unsigned long)> progressFunction;
		std::function<void(std::string&)> outputFunc;
		double estimatedVariance;
	};
	void generalisedSplittingFixedFactorsEvolution(generalisedSplittingFixedFactorsEvolutionArgs& args);
}
#endif
