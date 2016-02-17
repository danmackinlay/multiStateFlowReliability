#ifndef CRUDEMC_HEADER_GUARD
#define CRUDEMC_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
namespace multistateTurnip
{
	struct crudeMCArgs
	{
	public:
		crudeMCArgs(const Context& context)
			:context(context)
		{}
		const Context& context;
		double threshold;
		int n;
		int count;
		boost::mt19937 randomSource;
	};
	void crudeMC(crudeMCArgs& args);
}
#endif
