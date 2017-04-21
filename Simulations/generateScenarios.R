methods <- c("pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1")
#dodec5EqualCapacities <- expand.grid(method = c(methods, "gsFS"), demand = 5L, n = 50000L, nCapacities = 5, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18, 1e-19, 1e-20), graph = "dodecahedron5EqualCapacity")
dodec5EqualCapacities <- expand.grid(method = c(methods, "gsFS"), demand = 5L, n = 50000L, nCapacities = 5, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18), graph = "dodecahedron5EqualCapacity")

dodec15UnequalCapacities <- expand.grid(method = c(methods, "gsFS"), demand = 25L, n = 50000L, nCapacities = 15, interestVertices = "1,20", stringsAsFactors=FALSE, epsilon = c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13), graph = "dodecahedron15UnequalCapacity")

grid4x4_1Scenarios <- expand.grid(method = c("pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1", "gsFS"), demand = 20, n = 50000L, nCapacities = 8, interestVertices = "1,16", stringsAsFactors=FALSE, epsilon = c(1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13), graph = "grid4x4_1")

grid3x3_1Scenarios <- expand.grid(method = c("pmc", "turnipSingle", "turnipFull3", "turnipFull2", "turnipFull1", "gsFS"), demand = 20, n = 50000L, nCapacities = 8, interestVertices = "1,9", stringsAsFactors=FALSE, epsilon = c(1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13), graph = "grid3x3_1")

scenarios <- rbind(grid4x4_1Scenarios, grid3x3_1Scenarios, dodec15UnequalCapacities, dodec5EqualCapacities)

scenarios$file <- apply(scenarios, 1, function(x) paste0(as.numeric(x["epsilon"]), "-", as.integer(x["demand"]), "-", as.integer(x["n"]), "-", as.integer(x["nCapacities"]), "-", x["graph"], "-", x["method"], ".RData", sep=""))
scenarios$nReps <- 300
