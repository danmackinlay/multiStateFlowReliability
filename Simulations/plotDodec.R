library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tikzDevice)
colours <- brewer.pal(n = 7, name = "Dark2")
powerFormatter <- function(x)
{
	xAsCharacter <- as.character(x)
	#xAsCharacter <- sapply(xAsCharacter, function(y) if(nchar(y) == 2) paste0(y, ".0") else y)
	paste0("$10^{", xAsCharacter, "}$")
}
source("./generateScenarios.R")
load("./summarised.RData")

currentGraph <- "dodecahedron5EqualCapacity"
dir.create("dodecahedron5EqualCapacityPaper", showWarnings = FALSE)
setwd("dodecahedron5EqualCapacityPaper")
#Ignore the case of 0.1. I think there might be a bug with some of the methods in that case.....or the variances may be legitimately much lower, which is possible. 
currentGraphScenariosIndices <- which(scenarios$graph == currentGraph & scenarios$epsilon < 0.1)
currentGraphScenarios <- subset(scenarios, graph == currentGraph & epsilon < 0.1)
currentGraphScenarios$epsilon <- factor(as.character(currentGraphScenarios$epsilon), levels = as.character(sort(unique(currentGraphScenarios$epsilon), decreasing = TRUE)))
secondsPerRunData <- melt(cbind(currentGraphScenarios, secondsPerRun = averageSecondsPerRun[currentGraphScenariosIndices]), measure.vars = "secondsPerRun")
averageEstimatesData <- melt(cbind(currentGraphScenarios, averageEstimate = averageEstimates[currentGraphScenariosIndices]), measure.vars = "averageEstimate")
workNormalizedRelativeVariancePlusPilotData <- melt(cbind(currentGraphScenarios, workNormalizedRelativeVariance = workNormalizedRelativeVariancePlusPilot[currentGraphScenariosIndices]), measure.vars = "workNormalizedRelativeVariance")
nCapacityLevels <- unique(currentGraphScenarios$nCapacities)

for(nCapacityLevels in nCapacityLevels)
{
	demands <- unique(subset(currentGraphScenarios, nCapacities == nCapacityLevels)$demand)
	for(currentDemand in demands)
	{
		currentSecondsPerRunData <- subset(secondsPerRunData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("gsFS", "pmc", "turnipFull1", "turnipFull2"))
		timePlot <- ggplot(data=currentSecondsPerRunData, aes(x = epsilon, y = value, group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + ggtitle(paste0("Average time per 100000, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + ylab("Seconds") + xlab("epsilon") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9))
		pdf(paste0("./times-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(timePlot)
		dev.off()

		currentAverageEstimatesData <- subset(averageEstimatesData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("gsFS", "pmc", "turnipFull1", "turnipFull2"))
		estimatesPlot <- ggplot(data=currentAverageEstimatesData, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Estimate") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Average estimates, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels))
		pdf(paste0("./estimates-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(estimatesPlot)
		dev.off()
		
		currentWorkNormalizedRelativeVariancePlusPilotData <- subset(workNormalizedRelativeVariancePlusPilotData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("gsFS", "pmc", "turnipFull1", "turnipFull2"))
		workNormalizedRelativeVariancePlot <- ggplot(data=currentWorkNormalizedRelativeVariancePlusPilotData, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Work Normalized Relative Variance") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.25, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Work Normalized Relative Variances,  Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + theme(plot.title = element_text(size = rel(1)))
		pdf(paste0("./workNormalizedRelativeVariance-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(workNormalizedRelativeVariancePlot)
		dev.off()
		tikz(paste0("./workNormalizedRelativeVariance-", nCapacityLevels, "-", currentDemand, ".tex"))
			print(workNormalizedRelativeVariancePlot)	
		dev.off()
	}
}
setwd("..")
