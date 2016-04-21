pmc <- function(graph, capacityMatrix, n, threshold, seed, interestVertices, verbose=FALSE)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(vcount(graph), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("pmc_igraph", graph, capacityMatrix, n, threshold, seed, interestVertices, verbose, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("pmc_graphNEL", graph, capacityMatrix, n, threshold, seed, interestVertices, verbose, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("pmc_graphAM", graph, capacityMatrix, n, threshold, seed, interestVertices, verbose, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("pmcResult", call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, threshold = threshold, interestVertices = as.integer(interestVertices), seed = as.integer(seed), firstMomentSingleSample = mpfr(result$firstMomentSingleSample), secondMomentSingleSample = mpfr(result$secondMomentSingleSample), varianceSingleSample = mpfr(result$varianceSingleSample), sqrtVarianceOfEstimate = mpfr(result$sqrtVarianceOfEstimate), relativeErrorEstimate = mpfr(result$relativeErrorEstimate)))
}
