#ifndef CREATE_CAPACITY_DISTRIBUTION_HEADER_GUARD
#define CREATE_CAPACITY_DISTRIBUTION_HEADER_GUARD
#include "capacityDistribution.h"
#include "Rcpp.h"
namespace multistateTurnip
{
	capacityDistribution createCapacityDistribution(SEXP capacityData);
}
#endif
