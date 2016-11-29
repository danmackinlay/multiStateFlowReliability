source("./generateScenarios.R")
SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))


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
nReps <- scenarios[SCENARIO_INDEX, "nReps"]

cat("SCENARIO_INDEX=", SCENARIO_INDEX, "\n", sep="")
cat("method=", method, "\n", sep="")

getCapacityMatrix <- function(rho, bi, epsilon)
{
	capacityMatrix <- data.frame(capacity = 0:bi, probability = rho^(bi - (0:bi) - 1)*epsilon)
	capacityMatrix[bi+1,"probability"] <- 1 - sum(capacityMatrix[1:(bi),"probability"])
	return(capacityMatrix)
}

interestVertices <- as.integer(str_split(scenarios[SCENARIO_INDEX, "interestVertices"], ",")[[1]])

if(graph == "dodecahedron5EqualCapacity")
{
	graph <- igraph::read.graph("./dodecahedron.graphml", format = "graphml")
	capacityMatrix <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 4)
	capacityList <- replicate(30, capacityMatrix, simplify=FALSE)

	maxCapacity <- max(capacityMatrix[,1])
	maxPossibleFlow <- 3*maxCapacity
} else if(graph == "dodecahedron15UnequalCapacity")
{
	graph <- igraph::read.graph("./dodecahedron.graphml", format = "graphml")
	capacityMatrix1 <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 10)
	capacityMatrix2 <- getCapacityMatrix(rho = 0.7, epsilon = epsilon, bi = 15)
	capacityList <- replicate(30, capacityMatrix1, simplify=FALSE)

	maxPossibleFlow <- 30

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
	maxPossibleFlow <- 24
} else if(graph == "grid10x10_2")
{
	graph <- igraph::make_lattice(dimvector = c(10,10))
	graph <- igraph::add_vertices(graph, 2)
	extraEdges <- c(unlist(rbind(1:10, 101)), unlist(rbind(91:100, 102)))
	graph <- igraph::add_edges(graph, extraEdges)

	capacityMatrix1 <- getCapacityMatrix(rho = 0.5, epsilon = epsilon, bi = 10)
	capacityList <- replicate(180+20, capacityMatrix1, simplify=FALSE)
	maxPossibleFlow <- 100
} else if(graph == "Jane")
{
	janeEdges <- c(
		10, 1, 
		10, 4,
		10, 6, 
		10, 8, 
		1, 2, 
		4, 1,
		4, 5,
		6, 4, 
		6, 5,
		6, 7, 
		8, 6, 
		8, 7, 
		8, 9,
		5, 11, 
		7, 5, 
		2, 3, 
		11, 12, 
		3, 11, 
		9, 5, 
		9, 11, 
		4, 2)
	graph <- igraph::graph(edges = janeEdges, directed = TRUE)
	probabilitiesMatrix <- cbind(c(0.1163, 0.1624, 0.2014, 0.0689, 0.1863, 0.2244, 0.2221, 0.1265, 0.2993, 0.3016, 0.2385, 0.3459, 0.3511, 0.0326, 0.0231, 0.0373, 0.0222, 0.0052, 0.3935, 0.0650, 0.1260), c(0.0616, 0.1224, 0.0900, 0.1155, 0.1366, 0.0214, 0.1334, 0.0762, 0.0343, 0.0813, 0.0785, 0.0269, 0.0441, 0.0182, 0.1268, 0.0830, 0.0192, 0.0411, 0.0625, 0.0457, 0.0495), c(0.8221, 0.7152, 0.7086, 0.8156, 0.6771, 0.7542, 0.6445, 0.7973, 0.6664, 0.6171, 0.6830, 0.6272, 0.6048, 0.9492, 0.8501, 0.8797, 0.9586, 0.9537, 0.5440, 0.8893, 0.8245))
	capacityList <- lapply(1:nrow(probabilitiesMatrix), function(i) data.frame(capacity = c(0, 3, 5), probability = probabilitiesMatrix[i,]))
} else if(graph == "Alexopolous")
{
	edges <- c(
		1, 2, 
		1, 3,
		1, 4,
		2, 1, 
		2, 5, 
		2, 6, 
		3, 1, 
		3, 2, 
		3, 7, 
		3, 8, 
		4, 3, 
		4, 9,
		5, 7, 
		5, 10, 
		6, 3, 
		6, 5, 
		6, 7, 
		7, 6, 
		7, 8, 
		7, 10, 
		8, 4, 
		8, 7, 
		8, 9, 
		9, 7, 
		9, 10)
	graph <- igraph::graph(edges = edges, directed = TRUE)
	probabilitiesFunction <- function(x)
	{
		return(data.frame(capacity = c(0L, as.integer(x)), probability = c(0.1, 0.9)))
	}
	probabilitiesData <- lapply(as.list(c(10, 35, 24, 48, 50, 10, 45, 30, 24, 10, 12, 45, 18, 25, 12, 36, 50, 40, 15, 20, 45, 35, 25, 24, 30)), probabilitiesFunction)
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
	while(counter < nReps + 1)
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
	while(counter < nReps + 1)
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
	while(counter < nReps + 1)
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
	while(counter < nReps + 1)
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
	while(counter < nReps + 1)
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
	while(counter < nReps + 1)
	{
		results[[counter]] <- turnip(graph = graph, capacityMatrix = capacityList, n = n, threshold = demand, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, useAllPointsMaxFlow = TRUE, allPointsMaxFlowIncrement = 1L)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "gsFS")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else 
	{
		results <- list()
		pilot <- generalisedSplittingAdaptiveEvolution(graph = graph, capacityMatrix = capacityList, n = n, seed = SCENARIO_INDEX + counter * 100000 - 1, interestVertices = interestVertices, verbose=FALSE, fraction = 10, level = demand)
		factors <- rep(10, length(pilot@times)-1)
		save(pilot, factors, results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
	}
	while(counter < nReps + 1)
	{
		results[[counter]] <- generalisedSplittingFixedFactorsEvolution(graph = graph, capacityMatrix = capacityList, n = n, times = pilot@times, seed = SCENARIO_INDEX + counter * 100000, interestVertices = interestVertices, verbose=FALSE, factors = factors, level = demand)
		save(pilot, factors, results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else
{
	stop("Unrecognized method")
}
