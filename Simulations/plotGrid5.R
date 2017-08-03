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
powerFormatterLog <- function(x)
{
	xAsCharacter <- as.character(log10(as.numeric(x)))
	#xAsCharacter <- sapply(xAsCharacter, function(y) if(nchar(y) == 2) paste0(y, ".0") else y)
	paste0("$10^{", xAsCharacter, "}$")
}
source("./generateScenarios.R")
load("./summarised.RData")

currentGraph <- "grid5x5_1"
dir.create("grid5x5_1", showWarnings = FALSE)
setwd("grid5x5_1")

#Rename valus in method column
scenarios$method[scenarios$method == "gsFS"] <- "GS"
scenarios$method[scenarios$method == "pmc"] <- "PMC"
scenarios$method[scenarios$method == "turnipSingle"] <- "PMC-Single"
scenarios$method[scenarios$method == "turnipFull1"] <- "PMC-All"

currentGraphScenariosIndices <- which(scenarios$graph == currentGraph)
currentGraphScenarios <- subset(scenarios, graph == currentGraph & epsilon < 0.1)
currentGraphScenarios$epsilon <- factor(as.character(currentGraphScenarios$epsilon), levels = as.character(sort(unique(currentGraphScenarios$epsilon), decreasing = TRUE)))
secondsPerRunData <- melt(cbind(currentGraphScenarios, secondsPerRun = averageSecondsPerRun[currentGraphScenariosIndices]), measure.vars = "secondsPerRun")
averageEstimatesData <- melt(cbind(currentGraphScenarios, averageEstimate = averageEstimates[currentGraphScenariosIndices]), measure.vars = "averageEstimate")
workNormalizedRelativeVariancePlusPilotData <- melt(cbind(currentGraphScenarios, workNormalizedRelativeVariance = workNormalizedRelativeVariancePlusPilot[currentGraphScenariosIndices]), measure.vars = "workNormalizedRelativeVariance")
relativeErrorData <- melt(cbind(currentGraphScenarios, relativeError = relativeErrors[currentGraphScenariosIndices]), measure.vars = "relativeError")
nCapacityLevels <- unique(currentGraphScenarios$nCapacities)

for(nCapacityLevels in nCapacityLevels)
{
	demands <- unique(subset(currentGraphScenarios, nCapacities == nCapacityLevels)$demand)
	for(currentDemand in demands)
	{
		currentSecondsPerRunData <- subset(secondsPerRunData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("GS", "PMC", "PMC-All", "PMC-Single"))
		timePlot <- ggplot(data=currentSecondsPerRunData, aes(x = epsilon, y = value, group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + ggtitle(paste0("Average time per 100000, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels)) + ylab("Seconds") + xlab("epsilon") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9))
		pdf(paste0("./times-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(timePlot)
		dev.off()

		size <- 2
		currentAverageEstimatesData <- subset(averageEstimatesData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("GS", "PMC", "PMC-All", "PMC-Single"))
		estimatesPlot <- ggplot(data=currentAverageEstimatesData, aes(x = epsilon, y = log10(value), group = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("epsilon") + ylab("Estimate") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "epsilon") + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter) + ggtitle(paste0("Average estimates, Demand = ", currentDemand, ", Capacity levels = ", nCapacityLevels))
		pdf(paste0("./estimates-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(estimatesPlot)
		dev.off()
		
		currentWorkNormalizedRelativeVariancePlusPilotData <- subset(workNormalizedRelativeVariancePlusPilotData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("GS", "PMC", "PMC-All", "PMC-Single"))
		xBreaks <- unique(currentWorkNormalizedRelativeVariancePlusPilotData$epsilon)[seq(2, length(unique(currentWorkNormalizedRelativeVariancePlusPilotData$epsilon)), by = 2)]
		workNormalizedRelativeVariancePlot <- ggplot(data=currentWorkNormalizedRelativeVariancePlusPilotData, aes(x = epsilon, y = log10(value), group = method, linetype = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("$\\epsilon$") + ylab("Work Normalized Relative Variance") + geom_line(size = size) + scale_colour_manual(values = colours) + scale_x_discrete(name = "$\\epsilon$", labels = powerFormatterLog, breaks = xBreaks) + theme(legend.position = c(0.2, 0.8)) + scale_y_continuous(labels = powerFormatter) + ggtitle("") + theme(plot.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.25)), legend.key.width = unit(3, "line")) + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method"))
		pdf(paste0("./workNormalizedRelativeVariance-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(workNormalizedRelativeVariancePlot)
		dev.off()
		tikzFileName <- paste0("./workNormalizedRelativeVariance-", nCapacityLevels, "-", currentDemand, ".tex")
		tikz(tikzFileName, width = 7, height = 7, standAlone = TRUE, packages = c("\\usepackage{amssymb}", "\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}"))
			print(workNormalizedRelativeVariancePlot)	
		dev.off()
		system2("pdflatex", args = tikzFileName)

		currentRelativeErrorData <- subset(relativeErrorData, nCapacities == nCapacityLevels & demand == currentDemand & method %in% c("GS", "PMC", "PMC-All", "PMC-Single"))
		xBreaks <- unique(currentRelativeErrorData$epsilon)[seq(2, length(unique(currentRelativeErrorData$epsilon)), by = 2)]
		relativeErrorPlot <- ggplot(data=currentRelativeErrorData, aes(x = epsilon, y = log10(value), group = method, linetype = method, colour = method)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab("$\\epsilon$") + ylab("Relative Error") + geom_line(size = size) + scale_colour_manual(values = colours) + scale_x_discrete(name = "$\\epsilon$", labels = powerFormatterLog, breaks = xBreaks) + theme(legend.position = c(0.2, 0.85)) + scale_y_continuous(labels = powerFormatter) + ggtitle("") + theme(plot.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.25)), legend.key.width = unit(3, "line")) + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method"))
		pdf(paste0("./relativeError-", nCapacityLevels, "-", currentDemand, ".pdf"))
			print(relativeErrorPlot)
		dev.off()
		tikzFileName <- paste0("./relativeError-", nCapacityLevels, "-", currentDemand, ".tex")
		tikz(tikzFileName, width = 7, height = 7, standAlone = TRUE, packages = c("\\usepackage{amssymb}", "\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}"))
			print(relativeErrorPlot)	
		dev.off()
		system2("pdflatex", args = tikzFileName)
	}
}
setwd("..")
