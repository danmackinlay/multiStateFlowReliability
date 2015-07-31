#include "arguments.h"
#include <boost/iterator/counting_iterator.hpp>
#include <iostream>
namespace multistateTurnip
{
	bool readN(boost::program_options::variables_map& variableMap, int& out)
	{
		if(variableMap.count("n") != 1)
		{
			std::cout << "Please enter a single value for input `n'" << std::endl;
			return false;
		}
		out = variableMap["n"].as<int>();
		if(out <= 0)
		{
			std::cout << "Input `n' must be a positive integer" << std::endl;
			return false;
		}
		return true;
	}
	void readSeed(boost::program_options::variables_map& variableMap, boost::mt19937& randomSource)
	{
		if(variableMap.count("seed") > 0)
		{
			randomSource.seed(variableMap["seed"].as<int>());
		}
	}
}
