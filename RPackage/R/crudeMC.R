crudeMC <- function(graph, capacityMatrix, n, threshold, seed, interestVertices)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		start <- Sys.time()
		result <- .Call("crudeMC_igraph", graph, capacityMatrix, n, threshold, seed, interestVertices, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		start <- Sys.time()
		result <- .Call("crudeMC_graphNEL", graph, capacityMatrix, n, threshold, seed, interestVertices, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		start <- Sys.time()
		result <- .Call("crudeMC_graphAM", graph, capacityMatrix, n, threshold, seed, interestVertices, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("crudeMCResult", data = mpfr(result), call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, threshold = threshold, interestVertices = as.integer(interestVertices), seed = as.integer(seed)))
}
