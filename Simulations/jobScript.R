source("./generateScenarios.R")
SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
cat("SCENARIO_INDEX=", SCENARIO_INDEX, "\n", sep="")

nReps <- scenarios[SCENARIO_INDEX, "nReps"]

outputFile <- file.path("results", scenarios[SCENARIO_INDEX, "file"])
tmpFile <- paste0(outputFile, ".tmp")
library(multiStateFlowReliability)
library(stringr)
library(igraph)

epsilon <- scenarios[SCENARIO_INDEX, "epsilon"]
method <- scenarios[SCENARIO_INDEX, "method"]
demand <- scenarios[SCENARIO_INDEX, "demand"]
n <- scenarios[SCENARIO_INDEX, "n"]
graph <- scenarios[SCENARIO_INDEX, "graph"]
nCapacities <- scenarios[SCENARIO_INDEX, "nCapacities"]

getCapacityMatrix <- function(rho, bi, epsilon)
{
	capacityMatrix <- data.frame(capacity = 0:bi, probability = rho^(bi - (0:bi) - 1)*epsilon)
	capacityMatrix[bi+1,"probability"] <- 1 - sum(capacityMatrix[1:(bi),"probability"])
	return(capacityMatrix)
}

interestVertices <- as.integer(str_split(scenarios[SCENARIO_INDEX, "interestVertices"], ",")[[1]])

if(graph == "dodecahedron")
{
	graph <- igraph::read.graph("./dodecahedron.graphml", format = "graphml")
	capacityMatrix <- data.frame(capacity = 0:(nCapacities-1), probability = epsilon)
	capacityMatrix[nrow(capacityMatrix), "probability"] <- 1 - (nCapacities-1)*epsilon
	capacityList <- replicate(30, capacityMatrix, simplify=FALSE)

	maxCapacity <- max(capacityMatrix[,1])
	maxPossibleFlow <- 3*maxCapacity
	replications <- 100
} else if(graph == "grid10")
{
	graph <- igraph::graph.lattice(length = 10, dim = 2)

	capacityMatrix <- data.frame(capacity = 0:(nCapacities-1), probability = epsilon)
	capacityMatrix[nrow(capacityMatrix), "probability"] <- 1 - (nCapacities-1)*epsilon
	capacityList <- replicate(180, capacityMatrix, simplify=FALSE)

	maxCapacity <- max(capacityMatrix[,1])
	maxPossibleFlow <- 2*maxCapacity
	replications <- 100
} else if(graph == "dodecahedron5EqualCapacity")
{
	graph <- igraph::read.graph("./dodecahedron.graphml", format = "graphml")
	capacityMatrix <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 4)
	capacityList <- replicate(30, capacityMatrix, simplify=FALSE)

	maxCapacity <- max(capacityMatrix[,1])
	maxPossibleFlow <- 3*maxCapacity
	replications <- 1
} else if(graph == "dodecahedron15UnequalCapacity")
{
	graph <- igraph::read.graph("./dodecahedron.graphml", format = "graphml")
	capacityMatrix1 <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 10)
	capacityMatrix2 <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 15)
	capacityList <- replicate(30, capacityMatrix1, simplify=FALSE)

	maxPossibleFlow <- 30
	replications <- 1

	#Work out which edges should have the second capacity distribution
	edgeMatrix <- igraph::get.edges(graph, igraph::E(graph))
	secondCapacityEdges <- which(edgeMatrix[,1] %in% c(1, 20) | edgeMatrix[,2] %in% c(1, 20))
	capacityList[secondCapacityEdges] <- replicate(length(secondCapacityEdges), capacityMatrix2, simplify=FALSE)
} else if(graph == "grid10x10_1")
{
	graph <- igraph::make_lattice(dimvector = c(10,10))
	capacityMatrix1 <- getCapacityMatrix(rho = 0.6, epsilon = epsilon, bi = 8)
	capacityMatrix2 <- getCapacityMatrix(rho = 0.6, epsilon = epsilon, bi = 12)
	capacityList <- replicate(180, capacityMatrix1, simplify=FALSE)

	edgeMatrix <- igraph::get.edges(graph, igraph::E(graph))
	secondCapacityEdges <- which(edgeMatrix[,1] %in% c(1, 100) | edgeMatrix[,2] %in% c(1, 100))
	capacityList[secondCapacityEdges] <- replicate(length(secondCapacityEdges), capacityMatrix2, simplify=FALSE)
	replications <- 1
	maxPossibleFlow <- 24
} else
{
	stop("Unknown graph")
}

if(method == "crudeMC")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < replications + 1)
	{
		results[[counter]] <- crudeMC(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices) 
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "pmc")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < replications + 1)
	{
		results[[counter]] <- pmc(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000L, interestVertices = interestVertices)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "turnipSingle")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < replications + 1)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000L, interestVertices = interestVertices, useAllPointsMaxFlow = FALSE)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "turnipFull3")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < replications + 1)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000L, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 3L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "turnipFull2")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < replications + 1)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 2L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if (method == "turnipFull1")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < replications + 1)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 1L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "gsFE")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	if(maxPossibleFlow - demand %% 2 == 0)
	{
		levels <- seq(maxPossibleFlow, demand, -2)
	} else levels <- seq(maxPossibleFlow-1, demand, -2)
	while(counter < replications + 1)
	{
		results[[counter]] <- generalisedSplittingFixedEffort(graph = graph, capacityMatrix = capacityList, n = n, levels = levels, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, verbose=FALSE)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "gsFS")
{
	counter <- 1
	if(maxPossibleFlow - demand %% 2 == 0)
	{
		levels <- seq(maxPossibleFlow, demand, -2)
	} else levels <- seq(maxPossibleFlow-1, demand, -2)
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else 
	{
		results <- list()
		pilot <- generalisedSplittingFixedEffort(graph = graph, capacityMatrix = capacityList, n = n, levels = levels, seed = SCENARIO_INDEX + counter * 100000 - 1, interestVertices = interestVertices, verbose=FALSE)
		factors <- round(tail(1/pilot@levelProbabilities, -1))
	}
	while(counter < replications + 1)
	{
		results[[counter]] <- generalisedSplittingFixedFactors(graph = graph, capacityMatrix = capacityList, n = n, levels = levels, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, verbose=FALSE, factors = factors)
		save(pilot, factors, results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else
{
	stop("Unrecognized method")
}
