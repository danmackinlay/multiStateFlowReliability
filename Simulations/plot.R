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

setwd("plots")
graphs <- unique(scenarios$graph)
for(currentGraph in graphs)
{
	dir.create(currentGraph, showWarnings = FALSE)
	setwd(currentGraph)
	#Ignore the case of 0.1. I think there might be a bug with some of the methods in that case.....or the variances may be legitimately much lower, which is possible. 
	currentGraphScenariosIndices <- which(scenarios$graph == currentGraph & scenarios$epsilon < 0.1)
	currentGraphScenarios <- subset(scenarios, graph == currentGraph & epsilon < 0.1)
	currentGraphScenarios$epsilon <- factor(as.character(currentGraphScenarios$epsilon), levels = as.character(sort(unique(currentGraphScenarios$epsilon), decreasing = TRUE)))
	secondsPerRunData <- melt(cbind(currentGraphScenarios, secondsPerRun = averageSecondsPerRun[currentGraphScenariosIndices]), measure.vars = "secondsPerRun")
	secondsSingleSampleData <- melt(cbind(currentGraphScenarios, secondsSingleSample = secondsSingleSample[currentGraphScenariosIndices]), measure.vars = "secondsSingleSample")
	averageEstimatesData <- melt(cbind(currentGraphScenarios, averageEstimate = averageEstimates[currentGraphScenariosIndices]), measure.vars = "averageEstimate")
	workNormalizedSingleSampleVarianceData <- melt(cbind(currentGraphScenarios, workNormalizedSingleSampleVariance = workNormalizedSingleSampleVariance[currentGraphScenariosIndices]), measure.vars = "workNormalizedSingleSampleVariance")
	varianceDataSingleSample <- melt(cbind(currentGraphScenarios, varianceSingleSample = varianceSingleSample[currentGraphScenariosIndices]), measure.vars = "varianceSingleSample")
	nCapacityLevels <- unique(currentGraphScenarios$nCapacities)

	for(nCapacityLevels in nCapacityLevels)
	{
		demands <- unique(subset(currentGraphScenarios, nCapacities == nCapacityLevels)$demand)
		for(currentDemand in demands)
		{
			currentSecondsPerRunData <- subset(secondsPerRunData, nCapacities == nCapacityLevels & demand == currentDemand & method != "crudeMC")
			timePlot <- ggplot(data=currentSecondsPerRunData, aes(x = epsilon, y = value, group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + ggtitle(paste0("Average time per 100000, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + ylab("Seconds") + xlab("epsilon") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9))
			pdf(paste0("./times-", nCapacityLevels, "-", currentDemand, ".pdf"))
				print(timePlot)
			dev.off()

			currentSecondsSingleSampleData <- subset(secondsSingleSampleData, nCapacities == nCapacityLevels & demand == currentDemand & method != "crudeMC")
			timePlot <- ggplot(data=currentSecondsSingleSampleData, aes(x = epsilon, y = value, group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + ggtitle(paste0("Average time per sample, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + ylab("Seconds") + xlab("epsilon") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9))
			pdf(paste0("./timeSingleSample", nCapacityLevels, "-", currentDemand, ".pdf"))
				print(timePlot)
			dev.off()

			currentAverageEstimatesData <- subset(averageEstimatesData, nCapacities == nCapacityLevels & demand == currentDemand & method != "crudeMC")
			estimatesPlot <- ggplot(data=currentAverageEstimatesData, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Estimate") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Average estimates, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels))
			pdf(paste0("./estimates-", nCapacityLevels, "-", currentDemand, ".pdf"))
				print(estimatesPlot)
			dev.off()
			
			currentWorkNormalizedSingleSampleVarianceData <- subset(workNormalizedSingleSampleVarianceData, nCapacities == nCapacityLevels & demand == currentDemand & method != "crudeMC")
			workNormalizedSingleSampleVariancePlot <- ggplot(data=currentWorkNormalizedSingleSampleVarianceData, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Work Normalized Variance Per Sample") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Work Normalized Variances Per Sample, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + theme(plot.title = element_text(size = rel(1)))
			pdf(paste0("./workNormalizedSingleSampleVariance-", nCapacityLevels, "-", currentDemand, ".pdf"))
				print(workNormalizedSingleSampleVariancePlot)
			dev.off()

			currentWorkNormalizedRelativeVarianceSingleSampleData <- subset(workNormalizedSingleSampleVarianceData, nCapacities == nCapacityLevels & demand == currentDemand & method != "crudeMC")
			currentWorkNormalizedRelativeVarianceSingleSampleData$value <- currentWorkNormalizedRelativeVarianceSingleSampleData$value / (currentAverageEstimatesData$value^2)
			workNormalizedSingleSampleRelativeVariancePlot <- ggplot(data=currentWorkNormalizedRelativeVarianceSingleSampleData, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Work Normalized RelativeVariance Per Sample") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Work Normalized Relative Variances Per Sample, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + theme(plot.title = element_text(size = rel(1)))
			pdf(paste0("./workNormalizedSingleSampleRelativeVariance-", nCapacityLevels, "-", currentDemand, ".pdf"))
				print(workNormalizedSingleSampleRelativeVariancePlot)
			dev.off()

			currentVarianceDataSingleSample <- subset(varianceDataSingleSample, nCapacities == nCapacityLevels & demand == currentDemand & method != "crudeMC")
			varianceDataSingleSamplePlot <- ggplot(data=currentVarianceDataSingleSample, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Variance Per Sample") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Variances Per Sample, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels))
			pdf(paste0("./varianceSingleSample-", nCapacityLevels, "-", currentDemand, ".pdf"))
				print(varianceDataSingleSamplePlot)
			dev.off()
		}
	}
	setwd("..")
}
setwd("..")
