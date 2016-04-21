#include "crudeMC.h"
#include "createContext.h"
#include "createCapacityDistributions.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP crudeMC(SEXP graph, SEXP distributions_sexp, SEXP n_sexp, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, R_GRAPH_TYPE type)
	{
	BEGIN_RCPP
		double threshold;
		try
		{
			threshold = Rcpp::as<double>(threshold_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input threshold must be a single number");
		}

		std::vector<int> interestVertices;
		try
		{
			interestVertices = Rcpp::as<std::vector<int> >(interestVertices_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input interestVertices must be a pair of numbers");
		}

		int n;
		try
		{
			n = Rcpp::as<int>(n_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input n must be an integer");
		}

		int seed;
		try
		{
			seed = Rcpp::as<int>(seed_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input seed must be an integer");
		}
		if(interestVertices.size() != 2) std::runtime_error("Input interestVertices must be a pair of numbers");

		std::vector<capacityDistribution> distributions = createCapacityDistributions(distributions_sexp);
		Context context = createContext(graph, distributions, interestVertices[0]-1, interestVertices[1]-1, threshold, type);

		crudeMCArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.threshold = threshold;
		crudeMC(args);

		mpfr_class count = args.count, n_mpfr = args.n;
		mpfr_class result = count / n_mpfr;
		return Rcpp::wrap(result.str());
	END_RCPP
	}
	SEXP crudeMC_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return crudeMC(graph, distributions, n, threshold_sexp, seed_sexp, interestVertices_sexp, IGRAPH);
	}
	SEXP crudeMC_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return crudeMC(graph, distributions, n, threshold_sexp, seed_sexp, interestVertices_sexp, GRAPHAM);
	}
	SEXP crudeMC_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return crudeMC(graph, distributions, n, threshold_sexp, seed_sexp, interestVertices_sexp, GRAPHNEL);
	}

}
