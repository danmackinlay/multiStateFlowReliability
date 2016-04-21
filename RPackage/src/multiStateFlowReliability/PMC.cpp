#include "PMC.h"
#include "createContext.h"
#include "createCapacityDistributions.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP pmc(SEXP graph, SEXP distribution_sexp, SEXP n_sexp, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, R_GRAPH_TYPE type)
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

		bool verbose;
		try
		{
			verbose = Rcpp::as<bool>(verbose_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input verbose must be a boolean");
		}
		if(interestVertices.size() != 2) std::runtime_error("Input interestVertices must be a pair of numbers");

		Rcpp::RObject barHandle;
		Rcpp::Function txtProgressBar("txtProgressBar"), setTxtProgressBar("setTxtProgressBar"), close("close");
		std::function<void(unsigned long, unsigned long)> progressFunction = [](unsigned long, unsigned long){};
		if(verbose)
		{
			barHandle = txtProgressBar(Rcpp::Named("style") = 3, Rcpp::Named("min") = 0, Rcpp::Named("max") = 1000, Rcpp::Named("initial") = 0);
			progressFunction = [barHandle, setTxtProgressBar](unsigned long done, unsigned long totalSteps)
			{
				try
				{
					setTxtProgressBar.topLevelExec(barHandle, (int)((double)(1000*done) / (double)totalSteps));
				}
				catch(...){}
			};
		}

		std::vector<capacityDistribution> distribution = createCapacityDistributions(distribution_sexp);
		Context context = createContext(graph, distribution, interestVertices[0]-1, interestVertices[1]-1, threshold, type);

		pmcArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.threshold = threshold;
		args.outputFunc = [](std::string& output){Rcpp::Rcout << output << std::endl;};
		args.progressFunction = progressFunction;
		pmc(args);

		if(verbose)
		{
			close(barHandle);
		}

		std::string firstMomentSingleSampleStr = args.firstMomentSingleSample.str(), secondMomentSingleSampleStr = args.secondMomentSingleSample.str(), varianceSingleSampleStr = args.varianceSingleSample.str(), sqrtVarianceOfEstimateStr = args.sqrtVarianceOfEstimate.str(), relativeErrorEstimateStr = args.relativeErrorEstimate.str();
		return Rcpp::List::create(Rcpp::Named("firstMomentSingleSample") = firstMomentSingleSampleStr, Rcpp::Named("secondMomentSingleSample") = secondMomentSingleSampleStr, Rcpp::Named("varianceSingleSample") = varianceSingleSampleStr, Rcpp::Named("sqrtVarianceOfEstimate") = sqrtVarianceOfEstimateStr, Rcpp::Named("relativeErrorEstimate") = relativeErrorEstimateStr);
	END_RCPP
	}
	SEXP pmc_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp)
	{
		return pmc(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, IGRAPH);
	}
	SEXP pmc_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp)
	{
		return pmc(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, GRAPHAM);
	}
	SEXP pmc_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP threshold_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp)
	{
		return pmc(graph, capacity, n, threshold_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, GRAPHNEL);
	}

}
