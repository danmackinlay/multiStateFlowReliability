library(xtable)
source("./generateScenarios.R")
load("./summarised.RData")
graphs <- unique(scenarios$graph)
allResults <- cbind(scenarios, varianceSingleSample, workNormalizedSingleSampleVariancePlusPilot, workNormalizedSingleSampleVariance, secondsSingleSamplePlusPilot, secondsSingleSample, averageEstimates) 
captions <- list()
captions[["grid10x10_1"]] <- "Simulation results for the $10 \\times 10$ grid graph. The vertices of interest were diagonally opposite corners of the graph. For edges incident to the vertices of interest, $b_i = 12$. For all other edges, $b_i = 8$. $\\rho = 0.6$ for all edges. The demand threshold was 20. "
captions[["dodecahedron5EqualCapacity"]] <- "Simulation results for the dodecahedron graph. The vertices of interest were 1 and 20. Every edge had $b_i = 4$ and $\\rho = 0.7$. The demand threshold was 5. "
captions[["dodecahedron15UnequalCapacity"]] <- "Simulation results for the dodecahedron graph. The vertices of interest were 1 and 20. Every edge had $\\rho = 0.7$. The demand threshold was 25. Every edge incident to vertex 1 or 20 had $b_i = 15$. Every other edge had $b_i = 10$. "
captions[["grid10x10_2"]] <- "Simulation results for the augmented $10 \\times 10$ grid graph. The vertices of interest were vertices 101 and 102. Vertex 101 was connected to every vertex on one side of the grid, and vertex 102 connected to every vertex on the opposite side. For all edges $\\rho = 0.5$ and $b_i = 10$. The demand threshold was 50."
sink(file="tables.tex")
cat("\\documentclass{article}\n\\usepackage[a4paper,margin=1in,landscape]{geometry}\n\\begin{document}", sep="")
cat("Var is the variance of a single sample using the given method. $\\mathrm{WNV}_1$ is the work normalized variance per single sample, with the time for the pilot run included. $\\mathrm{WNV}_1$ is the work normalized variance per single sample, with the time for the pilot run excluded. $\\mathrm{T}_1$ is the total run-time, including the pilot run. $\\mathrm{T}_2$ is the total run-time, excluding the pilot run. $\\widehat{\\ell}$ is the estimate. Pilot runs used 50000 samples, and runs of all estimators used 5000000 samples.")

translation <- c("gsFS" = "GS", "turnipSingle" = "FilterJumps1", "turnipFull1" = "FilterJumps2 every step", "turnipFull2" = "FilterJumps2 every second step", "turnipFull3" = "FilterJumps2 every third step", "pmc" = "pmc", "crudeMC" = "crudeMC")
for(graph in graphs)
{
	subsetted <- allResults[allResults$graph == graph,]
	subsetted <- subsetted[,c("method", "epsilon", "varianceSingleSample", "workNormalizedSingleSampleVariancePlusPilot", "workNormalizedSingleSampleVariance", "secondsSingleSamplePlusPilot", "secondsSingleSample", "averageEstimates")]
	colnames(subsetted) <- c("Method", "$\\epsilon$", "$\\mathrm{Var}$", "$\\mathrm{WNV}_1$", "$\\mathrm{WNV}_2$", "$\\mathrm{T}_1$", "$\\mathrm{T}_2$", "$\\widehat{\\ell}$")
	subsetted$Method <- translation[subsetted$Method]
	table <- xtable(subsetted, display = c("s", "s", "e", "g", "g", "g", "g", "g", "g"), caption = captions[[graph]])
	print(table, include.rownames=FALSE, math.style.exponents=TRUE, sanitize.colnames.function=identity)
}
cat("\\end{document}", sep="")
sink()
