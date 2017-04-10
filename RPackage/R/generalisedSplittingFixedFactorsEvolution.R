generalisedSplittingFixedFactorsEvolution <- function(graph, capacityMatrix, n, times, level, seed, interestVertices, verbose=FALSE, factors)
{
	if(any(diff(times) >= 0)) stop("Input times must be decreasing")
	if(length(level) > 1) stop("Input level must be a single value")
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
		result <- .Call("generalisedSplittingFixedFactorsEvolution_igraph", graph, capacityMatrix, n, times, level, seed, interestVertices, verbose, factors, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("generalisedSplittingFixedFactorsEvolution_graphNEL", graph, capacityMatrix, n, times, level, seed, interestVertices, verbose, factors, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("generalisedSplittingFixedFactorsEvolution_graphAM", graph, capacityMatrix, n, times, level, seed, interestVertices, verbose, factors, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("generalisedSplittingFixedFactorsEvolutionResult", call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, times = times, level = level, interestVertices = as.integer(interestVertices), seed = as.integer(seed), estimate = result$estimate, factors = factors, timeProbabilities = result$timeProbabilities))
}
