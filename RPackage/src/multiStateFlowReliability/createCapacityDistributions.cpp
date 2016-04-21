#include "createCapacityDistributions.h"
namespace multistateTurnip
{
	std::vector<capacityDistribution> createCapacityDistributions(SEXP capacityList_sexp)
	{
		Rcpp::List capacityList;
		try
		{
			capacityList = Rcpp::as<Rcpp::List>(capacityList_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input capacityList must be a List");
		}
		std::vector<capacityDistribution> distributions;
		for(Rcpp::List::iterator entry = capacityList.begin(); entry != capacityList.end(); entry++)
		{
			Rcpp::DataFrame capacityData;
			try
			{
				capacityData = Rcpp::as<Rcpp::DataFrame>(*entry);
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
			distributions.emplace_back(capacityDistribution(capacityDataVector));
		}
		return distributions;
	}
}
