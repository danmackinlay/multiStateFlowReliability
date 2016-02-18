turnip <- function(graph, capacityMatrix, n, threshold, seed, interestVertices, useAllPointsMaxFlow, allPointsMaxFlowIncrement)
{
	if(missing(allPointsMaxFlowIncrement)) allPointsMaxFlowIncrement <- NA
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		start <- Sys.time()
		result <- .Call("turnip_igraph", graph, capacityMatrix, n, threshold, seed, interestVertices, useAllPointsMaxFlow, allPointsMaxFlowIncrement, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		start <- Sys.time()
		result <- .Call("turnip_graphNEL", graph, capacityMatrix, n, threshold, seed, interestVertices, useAllPointsMaxFlow, allPointsMaxFlowIncrement, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		start <- Sys.time()
		result <- .Call("turnip_graphAM", graph, capacityMatrix, n, threshold, seed, interestVertices, useAllPointsMaxFlow, allPointsMaxFlowIncrement, PACKAGE="multiStateFlowReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("turnipResult", call = match.call(), start = start, end = end, n = as.integer(n), capacity = capacityMatrix, threshold = threshold, interestVertices = as.integer(interestVertices), useAllPointsMaxFlow = useAllPointsMaxFlow, allPointsMaxFlowIncrement = as.integer(allPointsMaxFlowIncrement), seed = as.integer(seed), estimateFirstMoment = mpfr(result$estimateFirstMoment), estimateSecondMoment = mpfr(result$estimateSecondMoment), varianceEstimate = mpfr(result$varianceEstimate), sqrtVarianceEstimate = mpfr(result$sqrtVarianceEstimate), relativeErrorEstimate = mpfr(result$relativeErrorEstimate)))
}