methods <- c("crudeMC", "pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1")
epsilonValues <- c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
scenariosDodec3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L, 5L), n = 100000L, nCapacities = 3L, graph = "dodecahedron", method = methods, interestVertices = "1,20", stringsAsFactors = FALSE)
scenariosDodec10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 24L, 25L, 26L, 27L, 28L, 29L, 30L), n = 100000L, nCapacities = 11L, graph = "dodecahedron", method = methods, interestVertices = "1,20", stringsAsFactors = FALSE)

scenariosGrid10.3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L), n = 100000L, nCapacities = 3L, graph = "grid10", method = methods, interestVertices = "1,100", stringsAsFactors = FALSE)
scenariosGrid10.10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 14L, 15L, 16L, 17L, 18L, 19L, 20L), n = 100000L, nCapacities = 11L, graph = "grid10", method = methods, interestVertices = "1,100", stringsAsFactors = FALSE)

gsScenariosDodec3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L, 5L), n = 100000L, nCapacities = 3L, graph = "dodecahedron", method = c("gsFS", "gsFE"), interestVertices = "1,20", stringsAsFactors = FALSE)
gsScenariosDodec10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 24L, 25L, 26L, 27L, 28L, 29L, 30L), n = 100000L, nCapacities = 11L, graph = "dodecahedron", method = c("gsFS", "gsFE"), interestVertices = "1,20", stringsAsFactors = FALSE)
gsScenariosGrid10.3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L), n = 100000L, nCapacities = 3L, graph = "grid10", method = c("gsFS", "gsFE"), interestVertices = "1,100", stringsAsFactors = FALSE)
gsScenariosGrid10.10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 14L, 15L, 16L, 17L, 18L, 19L, 20L), n = 100000L, nCapacities = 11L, graph = "grid10", method = c("gsFS", "gsFE"), interestVertices = "1,100", stringsAsFactors = FALSE)
gsScenarios <- rbind(gsScenariosDodec3Levels, gsScenariosDodec10Levels, gsScenariosGrid10.3Levels, gsScenariosGrid10.10Levels)

dodec5EqualCapacities <- expand.grid(method = c(methods, c("gsFS", "gsFE")), demand = 5L, n = 1000000L, nCapacities = 5, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(0.01, 0.001, 0.0001, 0.00001, 0.000001), graph = "dodecahedron5EqualCapacity")

dodec15UnequalCapacities <- expand.grid(method = c(methods, c("gsFS", "gsFE")), demand = 25L, n = 1000000L, nCapacities = 15, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(0.01, 0.001, 0.0001, 0.00001, 0.000001), graph = "dodecahedron15UnequalCapacity")

grid1Scenarios <- data.frame(method = c(methods, c("gsFS", "gsFE")), demand = 20, n = 1000000L, nCapacities = 8, interestVertices = "1,100", stringsAsFactors=FALSE, epsilon = 0.01, graph = "grid10x10_1")

scenarios <- rbind(grid1Scenarios, dodec15UnequalCapacities, dodec5EqualCapacities, scenariosDodec3Levels, scenariosDodec10Levels, scenariosGrid10.3Levels, scenariosGrid10.10Levels, gsScenarios)

scenarios$file <- apply(scenarios, 1, function(x) paste0(as.numeric(x["epsilon"]), "-", as.integer(x["demand"]), "-", as.integer(x["n"]), "-", as.integer(x["nCapacities"]), "-", x["graph"], "-", x["method"], ".RData", sep=""))
