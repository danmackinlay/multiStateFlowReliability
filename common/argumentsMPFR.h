#ifndef ARGUMENTS_MPFR_HEADER_GUARD
#define ARGUMENTS_MPFR_HEADER_GUARD
#include <boost/program_options.hpp>
#include "includeMPFR.h"
#include "context.h"
#include "capacityDistribution.h"
namespace multistateTurnip
{
	bool readCapacityDistribution(boost::program_options::variables_map& variableMap, capacityDistribution& out, mpfr_class& originalMinFlow, std::string& error);
	bool readContext(boost::program_options::variables_map& variableMap, Context& out, capacityDistribution&& distribution, const mpfr_class& threshold);
	bool readThreshold(boost::program_options::variables_map& variableMap, mpfr_class& out);
}
#endif
