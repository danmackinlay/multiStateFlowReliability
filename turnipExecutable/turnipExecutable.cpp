#include <boost/program_options.hpp>
#include "arguments.h"
#include "argumentsMPFR.h"
#include "turnip.h"
#include <iostream>
namespace multistateTurnip
{
	int main(int argc, char **argv)
	{
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
			("full", boost::program_options::value<int>()->implicit_value(1), "(int) Run all-points max flow once for every n added edges. Defaults to 1. ")
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
		Context context = Context::emptyContext();
		double newThresholdD;
		{
			//Don't use this distribution - It gets moved into the Context object
			capacityDistribution movedDistribution;
			std::string error;
			//The original minimum possible flow value
			mpfr_class originalMinFlow;
			//Read the capacity distribution argument. It actually changes the distribution so that the minimum is zero, and returns the input minimum in originalMinFlow
			if(!readCapacityDistribution(variableMap, movedDistribution, originalMinFlow, error))
			{
				std::cout << error << std::endl;
				return 0;
			}
			//Change the threshold to account for the change in the minimum flow
			mpfr_class newThreshold = threshold - originalMinFlow;
			newThresholdD = (double)newThreshold;
			//truncate the set of possible capacities at the new threshold, otherwise sorting that many edge repair times is a bottleneck
			movedDistribution = movedDistribution.truncateAtMax(newThreshold);
			if(newThreshold < 0)
			{
				std::cout << "Minimum flow for all links is higher than the threshold" << std::endl;
				return 0;
			}
			if(!readContext(variableMap, context, std::move(movedDistribution), newThreshold))
			{
				return 0;
			}
		}
		bool fullSpecified = variableMap.count("full") > 0;
		int fullAllPointsIncrement = variableMap["full"].as<int>();
		turnipArgs arguments(context);
		arguments.n = n;
		arguments.threshold = newThresholdD;
		arguments.useAllPointsMaxFlow = fullSpecified;
		arguments.allPointsMaxFlowIncrement = fullAllPointsIncrement;

		readSeed(variableMap, arguments.randomSource);

		turnip(arguments);
	
		std::cout << "Unreliability probability estimate was " << (double)arguments.firstMomentSingleSample << std::endl;
		std::cout << "Relative error was " << (double)arguments.relativeErrorEstimate << std::endl;
		return 0;
	}
}

int main(int argc, char **argv)
{
	return multistateTurnip::main(argc, argv);
}
