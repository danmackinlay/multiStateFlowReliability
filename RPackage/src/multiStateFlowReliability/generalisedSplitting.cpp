#include "generalisedSplitting.h"
#include "createContext.h"
#include "createCapacityDistribution.h"
#include "Rcpp.h"
#include "convertGraph.h"
namespace multistateTurnip
{
	SEXP generalisedSplitting(SEXP graph, SEXP capacity, SEXP n_sexp, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, R_GRAPH_TYPE type)
	{
	BEGIN_RCPP
		std::vector<double> levels;
		try
		{
			levels = Rcpp::as<std::vector<double> >(levels_sexp);
		}
		catch(...)
		{
			throw std::runtime_error("Input threshold must be a single number");
		}
		for(int i = 1; i < (int)levels.size(); i++)
		{
			if(levels[i] >= levels[i-1]) throw std::runtime_error("Levels for generalised splitting must be decreasing");
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

		capacityDistribution distribution = createCapacityDistribution(capacity);
		Context context = createContext(graph, std::move(distribution), interestVertices[0]-1, interestVertices[1]-1, levels.back(), type);

		generalisedSplittingArgs args(context);
		args.randomSource.seed(seed);
		args.n = n;
		args.levels = levels;
		args.outputFunc = [](std::string& output){Rcpp::Rcout << output << std::endl;};
		args.progressFunction = progressFunction;
		generalisedSplitting(args);

		if(verbose)
		{
			close(barHandle);
		}
		return Rcpp::wrap(args.estimate);
	END_RCPP
	}
	SEXP generalisedSplitting_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp)
	{
		return generalisedSplitting(graph, capacity, n, levels_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, IGRAPH);
	}
	SEXP generalisedSplitting_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp)
	{
		return generalisedSplitting(graph, capacity, n, levels_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, GRAPHAM);
	}
	SEXP generalisedSplitting_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp)
	{
		return generalisedSplitting(graph, capacity, n, levels_sexp, seed_sexp, interestVertices_sexp, verbose_sexp, GRAPHNEL);
	}

}
