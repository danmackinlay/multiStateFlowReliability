crudeMC <- function(graph, capacityMatrix, n, threshold, seed, interestVertices, undirected)
{
	if(class(graph) == "igraph")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(vcount(graph), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("crudeMC_igraph", graph, capacityMatrix, n, threshold, seed, interestVertices, undirected, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("crudeMC_graphNEL", graph, capacityMatrix, n, threshold, seed, interestVertices, undirected, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		if(!is.list(capacityMatrix))
		{
			capacityMatrix <- replicate(length(nodes(graph)), capacityMatrix, simplify=FALSE)
		}
		start <- Sys.time()
		result <- .Call("crudeMC_graphAM", graph, capacityMatrix, n, threshold, seed, interestVertices, undirected, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("crudeMCResult", data = mpfr(result), call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, threshold = threshold, interestVertices = as.integer(interestVertices), seed = as.integer(seed), undirected = undirected))
}
