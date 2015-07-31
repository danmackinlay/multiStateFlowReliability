#include <boost/program_options.hpp>
#include "arguments.h"
#include "argumentsMPFR.h"
#include "includeMPFR.h"
#include <iostream>
namespace multistateTurnip
{
	int main(int argc, char **argv)
	{
		mpfr_set_default_prec(1024);

		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and completeGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph and completeGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("capacityFile", boost::program_options::value<std::string>(), "(string) File containing capacity distribution data. Incompatible with graphFile and gridGraph")
			("uniformCapacity", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) Pair of integers a and b, so that capacities are uniform on the integer range [a, b]")
			("threshold", boost::program_options::value<std::string>(), "(float) The threshold flow between the source and sink")
			("n", boost::program_options::value<int>(), "(int) The number of simulations to perform. ")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("help", "Display this message");
		
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << 
				"This program estimates the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component?\n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}

		int n;
		if(!readN(variableMap, n))
		{
			return 0;
		}
		mpfr_class threshold;
		if(!readThreshold(variableMap, threshold))
		{
			std::cout << "Unable to parse threshold value" << std::endl;
			return 0;
		}
		double thresholdD = (double)threshold;

		capacityDistribution distribution;
		std::string error;
		mpfr_class originalMinFlow;
		if(!readCapacityDistribution(variableMap, distribution, originalMinFlow, error))
		{
			std::cout << error << std::endl;
			return 0;
		}
		mpfr_class newThreshold = threshold - originalMinFlow;
		if(newThreshold < 0)
		{
			std::cout << "Minimum flow for all links is higher than the threshold" << std::endl;
			return 0;
		}
		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, std::move(distribution), newThreshold))
		{
			return 0;
		}

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);
		
		std::vector<double>& capacityVector = context.getCapacityVector();
		std::size_t nEdges = context.getNEdges();
		const Context::internalDirectedGraph& graph = context.getDirectedGraph();
		const capacityDistribution& randomCapacityDistribution = context.getDistribution();

		int count = 0;
		for(int i = 0; i < n; i++)
		{
			for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
			{
				capacityVector[2*edgeCounter] = capacityVector[2*edgeCounter + 1] = randomCapacityDistribution(randomSource);
			}
			double flow = context.getMaxFlow(capacityVector);
			if(flow < thresholdD) count++;
		}
		std::cout << "Unreliability probability estimate was " << ((double)count / (double)n) << std::endl;
		return 0;
	}
}

int main(int argc, char **argv)
{
	return multistateTurnip::main(argc, argv);
}
