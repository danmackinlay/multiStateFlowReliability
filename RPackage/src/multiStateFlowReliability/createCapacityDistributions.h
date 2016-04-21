#ifndef CREATE_CAPACITY_DISTRIBUTIONS_HEADER_GUARD
#define CREATE_CAPACITY_DISTRIBUTIONS_HEADER_GUARD
#include "capacityDistribution.h"
#include "Rcpp.h"
namespace multistateTurnip
{
	std::vector<capacityDistribution> createCapacityDistributions(SEXP capacityData);
}
#endif
