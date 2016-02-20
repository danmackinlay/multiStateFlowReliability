source("./generateScenarios.R")
SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
cat("SCENARIO_INDEX=", SCENARIO_INDEX, "\n", sep="")

nReps <- scenarios[SCENARIO_INDEX, "nReps"]

outputFile <- file.path("results", scenarios[SCENARIO_INDEX, "file"])
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
} else
{
	stop("Unknown graph")
}

if(method == "crudeMC")
{
	result <- crudeMC(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices) 
} else if(method == "pmc")
{
	result <- pmc(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices)
} else if(method == "turnipSingle")
{
	result <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices, useAllPointsMaxFlow = FALSE)
} else if(method == "turnipFull3")
{
	result <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 3L)
} else if(method == "turnipFull2")
{
	result <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 2L)
} else if (method == "turnipFull1")
{
	result <- turnip(graph = graph, capacityMatrix = capacityMatrix, n = n, threshold = demand, seed = SCENARIO_INDEX, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 1L)
} else
{
	stop("Unrecognized method")
}
save(result, file = outputFile)
