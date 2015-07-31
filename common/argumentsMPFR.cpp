#include "argumentsMPFR.h"
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
namespace multistateTurnip
{
	bool pairSorter(std::pair<mpfr_class, mpfr_class>& first, std::pair<mpfr_class, mpfr_class>& second)
	{
		return first.first < second.first;
	}
	bool readCapacityDistribution(boost::program_options::variables_map& variableMap, capacityDistribution& out, mpfr_class& originalMinFlow, std::string& error)
	{
		if (variableMap.count("capacityFile") + variableMap.count("uniformCapacity")  != 1)
		{
			error = "Please enter exactly one of `capacityFile' and `uniformCapacity'";
			return false;
		}
		std::vector<std::pair<mpfr_class, mpfr_class>> distributionData;
		if (variableMap.count("capacityFile") == 1)
		{
			std::string capacityFile = variableMap["capacityFile"].as<std::string>();
			std::ifstream fileHandle(capacityFile, std::ios_base::in);

			std::string line;
			mpfr_class sum = 0;
			while (std::getline(fileHandle, line))
			{
				std::vector<std::string> parts;
				boost::split(parts, line, boost::is_any_of(","), boost::token_compress_on);
				if (parts.size() != 2)
				{
					error = "Wrong number of entries on line";
					return false;
				}
				boost::trim(parts[0]);
				boost::trim(parts[1]);

				mpfr_class value, probability;
				if (mpfr_set_str(value.backend().data(), parts[0].c_str(), 10, MPFR_RNDN) != 0)
				{
					error = "Unable to parse first number";
					return false;
				}
				if (mpfr_set_str(probability.backend().data(), parts[1].c_str(), 10, MPFR_RNDN) != 0)
				{
					error = "Unable to parse second number";
					return false;
				}
				if (value < 0)
				{
					error = "Cannot have a negative flow";
					return false;
				}
				if (probability < 0 || probability > 1)
				{
					error = "Probabilities must be between 0 and 1";
					return false;
				}
				sum += probability;
				std::pair<mpfr_class, mpfr_class> dataPair(std::move(value), std::move(probability));
				distributionData.push_back(std::move(dataPair));
			}
			mpfr_class difference = (1 - sum);
			mpfr_class absDifference = boost::multiprecision::abs(difference);
			if (absDifference > 1e-6)
			{
				error = "Probabilities must sum to 1, to within an accuracy of 1e-6";
				return false;
			}

			std::sort(distributionData.begin(), distributionData.end(), pairSorter);
			originalMinFlow = distributionData[0].first;
			//Subtract original min flow away from all the stuff in the distribution data.
			for (std::vector<std::pair<mpfr_class, mpfr_class> >::iterator i = distributionData.begin(); i != distributionData.end(); i++)
			{
				i->first -= originalMinFlow;
			}
			out = capacityDistribution(distributionData);
		}
		else
		{
			std::vector<int> range = variableMap["uniformCapacity"].as<std::vector<int> >();
			if (range.size() != 2)
			{
				error = "Please enter exactly two integers if input uniformCapacity is used";
				return false;
			}
			if (range[0] < 0 || range[1] < 0)
			{
				error = "Integers entered for input uniformCapacity must be non-negative";
				return false;
			}
			if (range[0] >= range[1])
			{
				error = "Integers entered for input uniformCapacity must be form a range with at least two values";
				return false;
			}
			mpfr_class prob = 1 / mpfr_class(range[1] - range[0] + 1);
			originalMinFlow = range[0];
			for (int i = 0; i <= range[1] - range[0]; i++)
			{
				distributionData.push_back(std::pair<mpfr_class, mpfr_class>(i, prob));
			}
			out = capacityDistribution(distributionData);
		}
		return true;
	}
	bool readContext(boost::program_options::variables_map& variableMap, Context& out, capacityDistribution&& distribution, const mpfr_class& threshold)
	{
		boost::shared_ptr<std::vector<int> > interestVertices;
		if(variableMap.count("interestVertices") != 1)
		{
			std::cout << "Please enter only one value for input `interestVertices'" << std::endl;
			return false;
		}
		{
			std::vector<int> tmp = variableMap["interestVertices"].as<std::vector<int> >();
			interestVertices = boost::shared_ptr<std::vector<int> >(new std::vector<int>());
			interestVertices->insert(interestVertices->begin(), tmp.begin(), tmp.end());
		}

		int minInterest = *std::min_element(interestVertices->begin(), interestVertices->end());
		int maxInterest = *std::max_element(interestVertices->begin(), interestVertices->end());
		if(minInterest < 0)
		{
			std::cout << "Input `interestVertices' cannot contain negative indices" << std::endl;
			return false;
		}
		if(variableMap.count("graphFile") + variableMap.count("gridGraph") + variableMap.count("completeGraph") != 1)
		{
			std::cout << "Please enter exactly one of `completeGrahp', `gridGraph', `graphFile' or `torusGraph'" << std::endl;
			return false;
		}
		else if(variableMap.count("graphFile") == 1)
		{
			if(interestVertices->size() != 2)
			{
				std::cout << "Input `interestVertices' must contain two vertex indices" << std::endl;
				return false;
			}

			bool successful = false;
			std::string message;
			try
			{
				out = Context::fromFile(variableMap["graphFile"].as<std::string>(), successful, interestVertices, message, std::move(distribution), threshold);
			}
			catch(std::runtime_error& err)
			{
				message = err.what();
				successful = false;
			}
			if(!successful)
			{
				std::cout << "Error reading graphml file. " << message << ". Exiting..." << std::endl;
				return false;
			}
		}
		else if(variableMap.count("gridGraph") == 1)
		{
			if(interestVertices->size() != 2)
			{
				std::cout << "Input `interestVertices' must contain two vertex indices" << std::endl;
				return false;
			}

			int gridDimension = variableMap["gridGraph"].as<int>();
			if(gridDimension <= 0)
			{
				std::cout << "Input `gridGraph' must be a positive integer" << std::endl;
				return false;
			}
			if(maxInterest >= gridDimension*gridDimension || minInterest < 0)
			{
				std::cout << "Input 'interestVertices' must contain numbers between 0 and (nVertices - 1) inclusive" << std::endl;
				return false;
			}
			out = Context::gridContext(gridDimension, interestVertices, std::move(distribution), threshold);
		}
		else if(variableMap.count("completeGraph") == 1)
		{
			int nVertices = variableMap["completeGraph"].as<int>();
			if(interestVertices->size() != 1)
			{
				std::cout << "For a complete graph, input interestVertices must contain the number of vertices of interest" << std::endl;
				return false;
			}
			if(minInterest != 2)
			{
				std::cout << "There must be two vertices of interest" << std::endl;
				return false;
			}
			out = Context::completeContext(nVertices, minInterest, std::move(distribution), threshold);
		}
		return true;
	}
	bool readThreshold(boost::program_options::variables_map& variableMap, mpfr_class& out)
	{
		if (variableMap.count("threshold") != 1)
		{
			std::cout << "Please enter a single value for input `threshold'" << std::endl;
			return false;
		}
		int retVal = mpfr_set_str(out.backend().data(), variableMap["threshold"].as<std::string>().c_str(), 10, MPFR_RNDN);
		return retVal == 0;
	}
}
