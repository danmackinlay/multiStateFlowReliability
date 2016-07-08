#include "generalisedSplittingAdaptiveEvolution.h"
#include "createContext.h"
#include "createCapacityDistributions.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP generalisedSplittingAdaptiveEvolution(SEXP graph, SEXP distributions_sexp, SEXP n_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP fraction_sexp, R_GRAPH_TYPE type)
	{
	BEGIN_RCPP
		int fraction;
		try
		{
			fraction = Rcpp::as<int>(fraction_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input fraction must be an integer");
		}
		
		double level;
		try
		{
			level = Rcpp::as<double>(level_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input level must be a number");
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
		Context context = createContext(graph, distributions, interestVertices[0]-1, interestVertices[1]-1, level, type);

		generalisedSplittingAdaptiveEvolutionArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.level = level;
		args.outputFunc = [](std::string& output){Rcpp::Rcout << output << std::endl;};
		args.fraction = fraction;
		generalisedSplittingAdaptiveEvolution(args);

		return Rcpp::wrap(args.times);
	END_RCPP
	}
	SEXP generalisedSplittingAdaptiveEvolution_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp)
	{
		return generalisedSplittingAdaptiveEvolution(graph, capacity, n, level_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, factors_sexp, IGRAPH);
	}
	SEXP generalisedSplittingAdaptiveEvolution_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp)
	{
		return generalisedSplittingAdaptiveEvolution(graph, capacity, n, level_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, factors_sexp, GRAPHAM);
	}
	SEXP generalisedSplittingAdaptiveEvolution_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp)
	{
		return generalisedSplittingAdaptiveEvolution(graph, capacity, n, level_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, factors_sexp, GRAPHNEL);
	}

}
