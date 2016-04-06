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
		//There are two options, capacityFile and uniformCapacity
		if (variableMap.count("capacityFile") + variableMap.count("uniformCapacity")  != 1)
		{
			error = "Please enter exactly one of `capacityFile' and `uniformCapacity'";
			return false;
		}
		std::vector<std::pair<mpfr_class, mpfr_class> > distributionData;
		if (variableMap.count("capacityFile") == 1)
		{
			std::string capacityFile = variableMap["capacityFile"].as<std::string>();
			std::ifstream fileHandle(capacityFile, std::ios_base::in);

			std::string line;
			mpfr_class sum = 0;
			while (std::getline(fileHandle, line))
			{
				std::vector<std::string> parts;
				//Each line should be two numbers, comma separated
				boost::split(parts, line, boost::is_any_of(","), boost::token_compress_on);
				if (parts.size() != 2)
				{
					error = "Wrong number of entries on line";
					return false;
				}
				boost::trim(parts[0]);
				boost::trim(parts[1]);

				mpfr_class value, probability;
				//First value is the capacity value
				if (mpfr_set_str(value.backend().data(), parts[0].c_str(), 10, MPFR_RNDN) != 0)
				{
					error = "Unable to parse first number";
					return false;
				}
				//Second value is the probability of that capacity value
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
			//Check that the probabilities sum to 1
			mpfr_class difference = (1 - sum);
			mpfr_class absDifference = boost::multiprecision::abs(difference);
			if (absDifference > 1e-6)
			{
				error = "Probabilities must sum to 1, to within an accuracy of 1e-6";
				return false;
			}

			//Sort the data for the distribution
			std::sort(distributionData.begin(), distributionData.end(), pairSorter);
			//Alter the distribution so that the minimum value is zero, and return the original minimum in an output parameter
			originalMinFlow = distributionData[0].first;
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
		if(variableMap.count("source") != 1)
		{
			std::cout << "Please enter only one value for input `source'" << std::endl;
			return false;
		}
		if(variableMap.count("sink") != 1)
		{
			std::cout << "Please enter only one value for input `sink'" << std::endl;
			return false;
		}
		if(variableMap.count("graphFile") + variableMap.count("gridGraph") + variableMap.count("completeGraph") != 1)
		{
			std::cout << "Please enter exactly one of `completeGraph', `gridGraph', `graphFile' or `torusGraph'" << std::endl;
			return false;
		}
		else if(variableMap.count("completeGraph") == 1)
		{
			int nVertices = variableMap["completeGraph"].as<int>();
			out = Context::completeContext(nVertices, std::move(distribution), threshold);
		}
		else
		{
			int source = variableMap["source"].as<int>();
			int sink = variableMap["sink"].as<int>();

			if(source < 0 || sink < 0)
			{
				std::cout << "Inputs `source' and `sink' cannot be negative" << std::endl;
				return false;
			}

			if(variableMap.count("graphFile") == 1)
			{
				bool successful = false;
				std::string message;
				try
				{
					out = Context::fromFile(variableMap["graphFile"].as<std::string>(), successful, source, sink, message, std::move(distribution), threshold);
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
				int gridDimension = variableMap["gridGraph"].as<int>();
				if(gridDimension <= 0)
				{
					std::cout << "Input `gridGraph' must be a positive integer" << std::endl;
					return false;
				}
				if(source >= gridDimension*gridDimension || sink >= gridDimension*gridDimension)
				{
					std::cout << "Inputs `source' and `sink' must be between 0 and (nVertices - 1) inclusive" << std::endl;
					return false;
				}
				out = Context::gridContext(gridDimension, source, sink, std::move(distribution), threshold);
			}
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
