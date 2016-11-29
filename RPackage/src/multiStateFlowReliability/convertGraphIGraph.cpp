#include "convertGraph.h"
#include "context.h"
namespace multistateTurnip
{
	void convertGraphIGraph(SEXP graph_sexp, Context::inputGraph& graphRef)
	{
		//Convert graph object
		Rcpp::List graph;
		try
		{
			graph = Rcpp::as<Rcpp::List>(graph_sexp);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Unable to convert input graph to a list");
		}
		int nVertices = Rcpp::as<int>(graph(0));
		Rcpp::IntegerVector edgesVertex1 = Rcpp::as<Rcpp::IntegerVector>(graph(2));
		Rcpp::IntegerVector edgesVertex2 = Rcpp::as<Rcpp::IntegerVector>(graph(3));
		if(edgesVertex1.size() != edgesVertex2.size())
		{
			throw std::runtime_error("Internal error");
		}

		//Construct graph
		graphRef = Context::inputGraph(nVertices);
		for(int i = 0; i < edgesVertex1.size(); i++)
		{
			boost::add_edge(edgesVertex1(i), edgesVertex2(i), graphRef);
		}
	}
	void convertGraphDirectedIGraph(SEXP graph_sexp, ContextDirected::inputGraph& graphRef)
	{
		//Convert graph object
		Rcpp::List graph;
		try
		{
			graph = Rcpp::as<Rcpp::List>(graph_sexp);
		}
		catch(Rcpp::not_compatible&)
		{
			throw std::runtime_error("Unable to convert input graph to a list");
		}
		int nVertices = Rcpp::as<int>(graph(0));
		Rcpp::IntegerVector edgesVertex1 = Rcpp::as<Rcpp::IntegerVector>(graph(2));
		Rcpp::IntegerVector edgesVertex2 = Rcpp::as<Rcpp::IntegerVector>(graph(3));
		if(edgesVertex1.size() != edgesVertex2.size())
		{
			throw std::runtime_error("Internal error");
		}

		//Construct graph
		graphRef = ContextDirected::inputGraph(nVertices);
		int counter = 0;
		for(int i = 0; i < edgesVertex1.size(); i++)
		{
			boost::add_edge(edgesVertex1(i), edgesVertex2(i), counter, graphRef);
			counter++;
		}
	}
}
