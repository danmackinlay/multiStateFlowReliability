generalisedSplittingFixedFactors <- function(graph, capacityMatrix, n, levels, seed, interestVertices, verbose=FALSE, factors)
{
	if(any(diff(levels) >= 0)) stop("Input levels must be decreasing")
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
		result <- .Call("generalisedSplittingFixedFactors_igraph", graph, capacityMatrix, n, levels, seed, interestVertices, verbose, factors, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("generalisedSplittingFixedFactors_graphNEL", graph, capacityMatrix, n, levels, seed, interestVertices, verbose, factors, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("generalisedSplittingFixedFactors_graphAM", graph, capacityMatrix, n, levels, seed, interestVertices, verbose, factors, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("generalisedSplittingFixedFactorsResult", call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, levels = levels, interestVertices = as.integer(interestVertices), seed = as.integer(seed), estimate = result$estimate, factors = factors, levelProbabilities = result$levelProbabilities, estimatedVariance = result$estimatedVariance))
}
