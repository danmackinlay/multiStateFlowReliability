generalisedSplittingAdaptiveEvolution <- function(graph, capacityMatrix, n, level, seed, interestVertices, fraction, verbose=FALSE)
{
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
		result <- .Call("generalisedSplittingAdaptiveEvolution_igraph", graph, capacityMatrix, n, level, seed, interestVertices, verbose, fraction, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("generalisedSplittingAdaptiveEvolution_graphNEL", graph, capacityMatrix, n, level, seed, interestVertices, verbose, fraction, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("generalisedSplittingAdaptiveEvolution_graphAM", graph, capacityMatrix, n, level, seed, interestVertices, verbose, fraction, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("generalisedSplittingAdaptiveEvolutionResult", call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, level = level, interestVertices = as.integer(interestVertices), seed = as.integer(seed), fraction = fraction, times = result))
}
