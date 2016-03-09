source("./generateScenarios.R")
SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
cat("SCENARIO_INDEX=", SCENARIO_INDEX, "\n", sep="")

nReps <- scenarios[SCENARIO_INDEX, "nReps"]

outputFile <- file.path("results", scenarios[SCENARIO_INDEX, "file"])
tmpFile <- paste0(outputFile, ".tmp")
library(multiStateFlowReliability)
library(stringr)

epsilon <- scenarios[SCENARIO_INDEX, "epsilon"]
method <- scenarios[SCENARIO_INDEX, "method"]
demand <- scenarios[SCENARIO_INDEX, "demand"]
n <- scenarios[SCENARIO_INDEX, "n"]
graph <- scenarios[SCENARIO_INDEX, "graph"]
nCapacities <- scenarios[SCENARIO_INDEX, "nCapacities"]

interestVertices <- as.integer(str_split(scenarios[SCENARIO_INDEX, "interestVertices"], ",")[[1]])

capacityMatrix <- data.frame(capacity = 0:(nCapacities-1), probability = epsilon)
capacityMatrix[nrow(capacityMatrix), "probability"] <- 1 - (nCapacities-1)*epsilon

if(graph == "dodecahedron")
{
	graph <- igraph::read.graph("./dodecahedron.graphml", format = "graphml")
} else if(graph == "grid10")
{
	library(igraph)
	graph <- graph.lattice(length = 10, dim = 2)
} else
{
	stop("Unknown graph")
}

if(method == "crudeMC")
{
	counter <- 1
	results <- list()
	while(counter < 100)
	{
		results[[counter]] <- crudeMC(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices) 
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "pmc")
{
	counter <- 1
	results <- list()
	while(counter < 100)
	{
		results[[counter]] <- pmc(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000L, interestVertices = interestVertices)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "turnipSingle")
{
	counter <- 1
	results <- list()
	while(counter < 100)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000L, interestVertices = interestVertices, useAllPointsMaxFlow = FALSE)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "turnipFull3")
{
	counter <- 1
	results <- list()
	while(counter < 100)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000L, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 3L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "turnipFull2")
{
	counter <- 1
	results <- list()
	while(counter < 100)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 2L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if (method == "turnipFull1")
{
	counter <- 1
	results <- list()
	while(counter < 100)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 1L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else
{
	stop("Unrecognized method")
}
