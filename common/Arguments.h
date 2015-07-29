#ifndef ARGUMENTS_HEADER_GUARD
#define ARGUMENTS_HEADER_GUARD
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
namespace multistateTurnip
{
	bool readN(boost::program_options::variables_map& map, int& out);
	void readSeed(boost::program_options::variables_map& variableMap, boost::mt19937& randomSource);
}
#endif