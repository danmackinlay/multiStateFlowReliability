crudeMC <- function(graph, capacityMatrix, n, threshold, seed, interestVertices)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		result <- .Call("crudeMC_igraph", graph, capacityMatrix, n, threshold, seed, interestVertices, PACKAGE="multiStateFlowReliability")
	}
	else if(class(graph) == "graphNEL")
	{
		result <- .Call("crudeMC_graphNEL", graph, capacityMatrix, n, threshold, seed, interestVertices, PACKAGE="multiStateFlowReliability")
	}
	else if(class(graph) == "graphAM")
	{
		result <- .Call("crudeMC_graphAM", graph, capacityMatrix, n, threshold, seed, interestVertices, PACKAGE="multiStateFlowReliability")
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(result)
}
