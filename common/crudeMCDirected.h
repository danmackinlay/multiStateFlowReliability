#ifndef CRUDEMC_DIRECTED_HEADER_GUARD
#define CRUDEMC_DIRECTED_HEADER_GUARD
#include "contextDirected.h"
#include <boost/random/mersenne_twister.hpp>
namespace multistateTurnip
{
	struct crudeMCDirectedArgs
	{
	public:
		crudeMCDirectedArgs(const ContextDirected& context)
			:context(context)
		{}
		const ContextDirected& context;
		double threshold;
		int n;
		int count;
		boost::mt19937 randomSource;
	};
	void crudeMCDirected(crudeMCDirectedArgs& args);
}
#endif
