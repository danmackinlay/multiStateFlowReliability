pmc <- function(graph, capacityMatrix, n, removeRedundant, threshold, seed, interestVertices)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		start <- Sys.time()
		result <- .Call("pmc_igraph", graph, capacityMatrix, n, threshold, seed, interestVertices, removeRedundant, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		start <- Sys.time()
		result <- .Call("pmc_graphNEL", graph, capacityMatrix, n, threshold, seed, interestVertices, removeRedundant, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		start <- Sys.time()
		result <- .Call("pmc_graphAM", graph, capacityMatrix, n, threshold, seed, interestVertices, removeRedundant, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("pmcResult", call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, threshold = threshold, interestVertices = as.integer(interestVertices), removeRedundant = removeRedundant, seed = as.integer(seed), estimateFirstMoment = mpfr(result$estimateFirstMoment), estimateSecondMoment = mpfr(result$estimateSecondMoment), varianceEstimate = mpfr(result$varianceEstimate), sqrtVarianceEstimate = mpfr(result$sqrtVarianceEstimate), relativeErrorEstimate = mpfr(result$relativeErrorEstimate)))
}