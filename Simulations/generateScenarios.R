methods <- c("crudeMC", "pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1")
dodec5EqualCapacities <- expand.grid(method = c(methods, "gsFS"), demand = 5L, n = 50000L, nCapacities = 5, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(0.01, 0.001, 0.0001, 1e-5, 1e-6), graph = "dodecahedron5EqualCapacity")

dodec15UnequalCapacities <- expand.grid(method = c(methods, "gsFS"), demand = 25L, n = 50000L, nCapacities = 15, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(0.01, 0.001, 0.0001, 1e-5, 1e-6), graph = "dodecahedron15UnequalCapacity")

grid10x10_1Scenarios <- expand.grid(method = c("pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1", "gsFS"), demand = 20, n = 50000L, nCapacities = 8, interestVertices = "1,100", stringsAsFactors=FALSE, epsilon = c(0.01, 0.001, 0.0001, 0.00001, 0.000001, 1e-7, 1e-8), graph = "grid10x10_1")

grid10x10_2Scenarios <- expand.grid(method = c("pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1", "gsFS"), demand = 50, n = 50000L, nCapacities = 10, interestVertices = "101,102", stringsAsFactors=FALSE, epsilon = 0.01, graph = "grid10x10_2")

janeScenarios <- expand.grid(method = c("pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1", "gsFS"), demand = 1:5, n = 50000L, nCapacities = -1, interestVertices = "10,12", stringsAsFactors = FALSE, epsilon = -1, graph = "Jane")

scenarios <- rbind(grid10x10_1Scenarios, grid10x10_2Scenarios, dodec15UnequalCapacities, dodec5EqualCapacities, janeScenarios)

scenarios$file <- apply(scenarios, 1, function(x) paste0(as.numeric(x["epsilon"]), "-", as.integer(x["demand"]), "-", as.integer(x["n"]), "-", as.integer(x["nCapacities"]), "-", x["graph"], "-", x["method"], ".RData", sep=""))
scenarios$nReps <- 100
