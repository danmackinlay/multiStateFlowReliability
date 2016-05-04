context("estimate variance of the gsFS estimator")

test_that("Check that the gsFS variance estimation is accurate",
{
	getCapacityMatrix <- function(rho, bi, epsilon)
	{
		capacityMatrix <- data.frame(capacity = 0:bi, probability = rho^(bi - (0:bi) - 1)*epsilon)
		capacityMatrix[bi+1,"probability"] <- 1 - sum(capacityMatrix[1:(bi),"probability"])
		return(capacityMatrix)
	}

	epsilon <- 0.01
	demand <- 25
	graphFile <- system.file("data", "dodecahedron.graphml", package="multiStateFlowReliability")
	graph <- igraph::read.graph(graphFile, format = "graphml")
	capacityMatrix <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 10)
	capacityList <- replicate(30, capacityMatrix, simplify=FALSE)

	maxPossibleFlow <- 30
	replications <- 400

	counter <- 1
	if(maxPossibleFlow - demand %% 2 == 0)
	{
		levels <- seq(maxPossibleFlow, demand, -2)
	} else levels <- seq(maxPossibleFlow-1, demand, -2)

	results <- list()
	pilot <- generalisedSplittingFixedEffort(graph = graph, capacityMatrix = capacityList, n = 50000, levels = levels, seed = 1, interestVertices = c(1, 20), verbose=FALSE)
		factors <- round(tail(1/pilot@levelProbabilities, -1))
		while(counter < replications + 1)
	{
		results[[counter]] <- generalisedSplittingFixedFactors(graph = graph, capacityMatrix = capacityList, n = 50000, levels = levels, seed = counter, interestVertices = c(1, 20), verbose=FALSE, factors = factors)
		counter <- counter + 1
	}
	averageOfVariances <- mean(unlist(lapply(results, function(x) x@estimatedVariance)))
	varianceOfEstimates <- var(unlist(lapply(results, function(x) x@estimate)))
	expect_equal(averageOfVariances, varianceOfEstimates, tolerance=0.05)
})
