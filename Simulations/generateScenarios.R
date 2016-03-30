methods <- c("crudeMC", "pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1")
epsilonValues <- c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
scenariosDodec3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L, 5L), n = 100000L, nCapacities = 3L, graph = "dodecahedron", method = methods, interestVertices = "1,20", stringsAsFactors = FALSE)
scenariosDodec10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 24L, 25L, 26L, 27L, 28L, 29L, 30L), n = 100000L, nCapacities = 11L, graph = "dodecahedron", method = methods, interestVertices = "1,20", stringsAsFactors = FALSE)

scenariosGrid10.3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L), n = 100000L, nCapacities = 3L, graph = "grid10", method = methods, interestVertices = "1,100", stringsAsFactors = FALSE)
scenariosGrid10.10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 14L, 15L, 16L, 17L, 18L, 19L, 20L), n = 100000L, nCapacities = 11L, graph = "grid10", method = methods, interestVertices = "1,100", stringsAsFactors = FALSE)

gsScenariosDodec3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L, 5L), n = 100000L, nCapacities = 3L, graph = "dodecahedron", method = "gs", interestVertices = "1,20", stringsAsFactors = FALSE)
gsScenariosDodec10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 24L, 25L, 26L, 27L, 28L, 29L, 30L), n = 100000L, nCapacities = 11L, graph = "dodecahedron", method = "gs", interestVertices = "1,20", stringsAsFactors = FALSE)
gsScenariosGrid10.3Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 4L), n = 100000L, nCapacities = 3L, graph = "grid10", method = "gs", interestVertices = "1,100", stringsAsFactors = FALSE)
gsScenariosGrid10.10Levels <- expand.grid(epsilon = epsilonValues, demand = c(1L, 2L, 3L, 14L, 15L, 16L, 17L, 18L, 19L, 20L), n = 100000L, nCapacities = 11L, graph = "grid10", method = "gs", interestVertices = "1,100", stringsAsFactors = FALSE)
gsScenarios <- rbind(gsScenariosDodec3Levels, gsScenariosDodec10Levels, gsScenariosGrid10.3Levels, gsScenariosGrid10.10Levels)

scenarios <- rbind(scenariosDodec3Levels, scenariosDodec10Levels, scenariosGrid10.3Levels, scenariosGrid10.10Levels, gsScenarios)

scenarios$file <- apply(scenarios, 1, function(x) paste0(as.numeric(x["epsilon"]), "-", as.integer(x["demand"]), "-", as.integer(x["n"]), "-", as.integer(x["nCapacities"]), "-", x["graph"], "-", x["method"], ".RData", sep=""))
