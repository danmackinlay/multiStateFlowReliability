#include "createCapacityDistribution.h"
namespace multistateTurnip
{
	capacityDistribution createCapacityDistribution(SEXP capacityData_sexp)
	{
		Rcpp::DataFrame capacityData;
		try
		{
			capacityData = Rcpp::as<Rcpp::DataFrame>(capacityData_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input capacityData must be a data.frame");
		}
		Rcpp::IntegerVector capacities;
		Rcpp::NumericVector probabilities;
		try
		{
			capacities = capacityData("capacity");
		}
		catch(...)
		{
			throw std::runtime_error("Column capacity must be an integer vector");
		}
		try
		{
			probabilities = capacityData("probability");
		}
		catch(...)
		{
			throw std::runtime_error("Column probability must be an integer vector");
		}

		std::vector<std::pair<mpfr_class, mpfr_class> > capacityDataVector;
		for(int i = 0; i < capacityData.nrows(); i++)
		{
			mpfr_class capacity = capacities(i);
			mpfr_class probability = probabilities(i);
			capacityDataVector.push_back(std::make_pair(capacity, probability));
		}
		return capacityDistribution(capacityDataVector);
	}
}
