#include "crudeMC.h"
namespace multistateTurnip
{
	void crudeMC(crudeMCArgs& args)
	{
		const Context& context = args.context;
		std::vector<double>& capacityVector = context.getCapacityVector();
		std::size_t nEdges = context.getNEdges();
		const capacityDistribution& randomCapacityDistribution = context.getDistribution();

		args.count = 0;
		for(int i = 0; i < args.n; i++)
		{
			for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
			{
				capacityVector[2*edgeCounter] = capacityVector[2*edgeCounter + 1] = randomCapacityDistribution(args.randomSource);
			}
			double flow = context.getMaxFlow(capacityVector);
			if(flow <= args.threshold) args.count++;
		}
	}
}
