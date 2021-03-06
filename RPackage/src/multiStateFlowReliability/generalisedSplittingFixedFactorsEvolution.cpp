#include "generalisedSplittingFixedFactorsEvolution.h"
#include "createContext.h"
#include "createCapacityDistributions.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP generalisedSplittingFixedFactorsEvolution(SEXP graph, SEXP distributions_sexp, SEXP n_sexp, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp, R_GRAPH_TYPE type)
	{
	BEGIN_RCPP
		std::vector<double> times;
		try
		{
			times = Rcpp::as<std::vector<double> >(times_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input times must be a numeric vector");
		}
		for(int i = 1; i < (int)times.size(); i++)
		{
			if(times[i] >= times[i-1]) throw std::runtime_error("Times for generalised splitting must be decreasing");
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

		std::vector<int> factors;
		try
		{
			factors = Rcpp::as<std::vector<int> >(factors_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input factors must be an integer vector");
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

		std::vector<capacityDistribution> distributions = createCapacityDistributions(distributions_sexp);
		Context context = createContext(graph, distributions, interestVertices[0]-1, interestVertices[1]-1, level, type);

		generalisedSplittingFixedFactorsEvolutionArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.times = times;
		args.level = level;
		args.outputFunc = [](std::string& output){Rcpp::Rcout << output << std::endl;};
		args.progressFunction = progressFunction;
		args.splittingFactors = factors;
		generalisedSplittingFixedFactorsEvolution(args);

		if(verbose)
		{
			close(barHandle);
		}
		return Rcpp::List::create(Rcpp::Named("estimate") = args.estimate, Rcpp::Named("timeProbabilities") = args.timeProbabilities);
	END_RCPP
	}
	SEXP generalisedSplittingFixedFactorsEvolution_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp)
	{
		return generalisedSplittingFixedFactorsEvolution(graph, capacity, n, times_sexp, level_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, factors_sexp, IGRAPH);
	}
	SEXP generalisedSplittingFixedFactorsEvolution_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp)
	{
		return generalisedSplittingFixedFactorsEvolution(graph, capacity, n, times_sexp, level_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, factors_sexp, GRAPHAM);
	}
	SEXP generalisedSplittingFixedFactorsEvolution_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp)
	{
		return generalisedSplittingFixedFactorsEvolution(graph, capacity, n, times_sexp, level_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, factors_sexp, GRAPHNEL);
	}

}
