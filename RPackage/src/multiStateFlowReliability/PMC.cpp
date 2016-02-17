#include "PMC.h"
#include "createContext.h"
#include "createCapacityDistribution.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP pmc(SEXP graph, SEXP capacity, SEXP n_sexp, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP removeRedundant_sexp, R_GRAPH_TYPE type)
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

		bool removeRedundant;
		try
		{
			removeRedundant = Rcpp::as<bool>(removeRedundant_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input removeRedundant must be a boolean");
		}

		capacityDistribution distribution = createCapacityDistribution(capacity);
		Context context = createContext(graph, std::move(distribution), interestVertices[0]-1, interestVertices[1]-1, threshold, type);

		pmcArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.threshold = threshold;
		args.useTurnip = removeRedundant;
		args.outputFunc = [](std::string& output){Rcpp::Rcout << output << std::endl;};
		pmc(args);

		std::string estimateFirstMomentStr = args.estimateFirstMoment.str(), estimateSecondMomentStr = args.estimateSecondMoment.str(), varianceEstimateStr = args.varianceEstimate.str(), sqrtVarianceEstimateStr = args.sqrtVarianceEstimate.str(), relativeErrorEstimateStr = args.relativeErrorEstimate.str();
		return Rcpp::List::create(Rcpp::Named("estimateFirstMoment") = estimateFirstMomentStr, Rcpp::Named("estimateSecondMoment") = estimateSecondMomentStr, Rcpp::Named("varianceEstimate") = varianceEstimateStr, Rcpp::Named("sqrtVarianceEstimate") = sqrtVarianceEstimateStr, Rcpp::Named("relativeErrorEstimate") = relativeErrorEstimateStr);
	END_RCPP
	}
	SEXP pmc_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP removeRedundant_sexp)
	{
		return pmc(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, removeRedundant_sexp, IGRAPH);
	}
	SEXP pmc_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP removeRedundant_sexp)
	{
		return pmc(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, removeRedundant_sexp, GRAPHAM);
	}
	SEXP pmc_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP removeRedundant_sexp)
	{
		return pmc(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, removeRedundant_sexp, GRAPHNEL);
	}

}
