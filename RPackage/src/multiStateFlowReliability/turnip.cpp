#include "turnip.h"
#include "createContext.h"
#include "createCapacityDistributions.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP turnip(SEXP graph, SEXP distributions_sexp, SEXP n_sexp, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp, R_GRAPH_TYPE type)
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

		bool useAllPointsMaxFlow;
		try
		{
			useAllPointsMaxFlow = Rcpp::as<bool>(useAllPointsMaxFlow_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input useAllPointsMaxFlow must be a boolean");
		}
		int allPointsMaxFlowIncrement = 0;
		if(useAllPointsMaxFlow)
		{
			try
			{
				allPointsMaxFlowIncrement = Rcpp::as<int>(allPointsMaxFlowIncrement_sexp);
			}
			catch(...)
			{
				throw std::runtime_error("Input allPointsMaxFlowIncrement must be an integer");
			}
		}


		std::vector<capacityDistribution> distributions = createCapacityDistributions(distributions_sexp);
		Context context = createContext(graph, distributions, interestVertices[0]-1, interestVertices[1]-1, threshold, type);

		turnipArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.threshold = threshold;
		args.useAllPointsMaxFlow = useAllPointsMaxFlow;
		args.allPointsMaxFlowIncrement = allPointsMaxFlowIncrement;
		turnip(args);

		std::string firstMomentSingleSampleStr = args.firstMomentSingleSample.str(), secondMomentSingleSampleStr = args.secondMomentSingleSample.str(), varianceSingleSampleStr = args.varianceSingleSample.str(), sqrtVarianceOfEstimateStr = args.sqrtVarianceOfEstimate.str(), relativeErrorEstimateStr = args.relativeErrorEstimate.str();
		return Rcpp::List::create(Rcpp::Named("firstMomentSingleSample") = firstMomentSingleSampleStr, Rcpp::Named("secondMomentSingleSample") = secondMomentSingleSampleStr, Rcpp::Named("varianceSingleSample") = varianceSingleSampleStr, Rcpp::Named("sqrtVarianceOfEstimate") = sqrtVarianceOfEstimateStr, Rcpp::Named("relativeErrorEstimate") = relativeErrorEstimateStr);
	END_RCPP
	}
	SEXP turnip_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp)
	{
		return turnip(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, useAllPointsMaxFlow_sexp, allPointsMaxFlowIncrement_sexp, IGRAPH);
	}
	SEXP turnip_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp)
	{
		return turnip(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, useAllPointsMaxFlow_sexp, allPointsMaxFlowIncrement_sexp, GRAPHAM);
	}
	SEXP turnip_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP useAllPointsMaxFlow_sexp, SEXP allPointsMaxFlowIncrement_sexp)
	{
		return turnip(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, useAllPointsMaxFlow_sexp, allPointsMaxFlowIncrement_sexp, GRAPHNEL);
	}

}
