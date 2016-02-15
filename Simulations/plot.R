load("./summarised.RData")
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

scenarios$Distribution <- as.factor(scenarios$Distribution)
timeData <- melt(cbind(scenarios, times), measure.vars = colnames(times))
timePerSampleData <- melt(cbind(scenarios, timesPerSample), measure.vars = colnames(times))
estimatesData <- melt(cbind(scenarios, estimates), measure.vars = colnames(estimates))
relativeErrorData <- melt(cbind(scenarios, relativeErrors), measure.vars = colnames(relativeErrors))

workNormalizedRelativeErrorData <- relativeErrorData
workNormalizedRelativeErrorData$value <- workNormalizedRelativeErrorData$value * sqrt(timeData$value)

levels(relativeErrorData$variable) <- levels(estimatesData$variable) <- levels(timeData$variable) <- colnames(times)

thresholdValues <- unique(scenarios[,"Threshold"])
distributionValues <- unique(scenarios[,"Distribution"])
for(threshold in unique(thresholdValues))
{
	currentThresholdTimeData <- subset(timeData, Threshold == threshold)
	timePlot <- ggplot(data=currentThresholdTimeData, aes(x = Distribution, y = value, group = variable, colour = variable)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab(paste0("Threshold = ", threshold)) + ylab("Seconds") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "Distribution", labels = distributionValues) + theme(legend.position = c(0.77, 0.9))
	pdf(paste0("./times", threshold, ".pdf"))
		print(timePlot)
	dev.off()

	currentThresholdEstimatesData <- subset(estimatesData, Threshold == threshold)
	estimatesPlot <- ggplot(data=currentThresholdEstimatesData, aes(x = Distribution, y = log10(value), group = variable, colour = variable)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab(paste0("Threshold = ", threshold)) + ylab("Estimate") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "Distribution", labels = distributionValues) + theme(legend.position = c(0.77, 0.9)) + scale_y_continuous(labels = powerFormatter)
	texFile <- paste0("./estimates", threshold, ".tex")
	tikz(texFile, width = 7, height = 7, standAlone = TRUE, packages = c("\\usepackage{amssymb}", "\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}"))
		print(estimatesPlot)
	dev.off()
	system(paste0("pdflatex ", texFile))
	unlink(paste0("./estimates", threshold, ".tex"))
	unlink(paste0("./estimates", threshold, ".aux"))
	unlink(paste0("./estimates", threshold, ".log"))

	currentThresholdTimePerSampleData <- subset(timePerSampleData, Threshold == threshold)
	timesPerSamplePlot <- ggplot(data=currentThresholdTimePerSampleData, aes(x = Distribution, y = value, group = variable, colour = variable)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab(paste0("Threshold = ", threshold)) + ylab("Seconds per sample") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "Distribution", labels = distributionValues) + theme(legend.position = c(0.77, 0.9))
	pdf(paste0("./timesPerSample", threshold, ".pdf"))
		print(timesPerSamplePlot)
	dev.off()

	currentThresholdRelativeError <- subset(relativeErrorData, Threshold == threshold)
	relativeErrorPlot <- ggplot(data=currentThresholdRelativeError, aes(x = Distribution, y = value, group = variable, colour = variable)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab(paste0("Threshold = ", threshold)) + ylab("Relative error") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "Distribution", labels = distributionValues) + theme(legend.position = c(0.23, 0.895)) + coord_cartesian(ylim = c(0, 1))
	pdf(paste0("./relativeError", threshold, ".pdf"))
		print(relativeErrorPlot)
	dev.off()

	currentThresholdWorkNormalizedRelativeError <- subset(workNormalizedRelativeErrorData, Threshold == threshold)
	workNormalizedRelativeErrorPlot <- ggplot(data=currentThresholdWorkNormalizedRelativeError, aes(x = Distribution, y = value, group = variable, colour = variable)) + guides(col=guide_legend(ncol = 2, title=NULL)) + xlab(paste0("Threshold = ", threshold)) + ylab("Work normalized relative error") + geom_line() + scale_colour_manual(values = colours) + scale_x_discrete(name = "Distribution", labels = distributionValues) + theme(legend.position = c(0.23, 0.895)) + coord_cartesian(ylim = c(0, 1))
	pdf(paste0("./workNormalizedRelativeError", threshold, ".pdf"))
		print(workNormalizedRelativeErrorPlot)
	dev.off()
}

estimatesXTable <- print(xtable(cbind(scenarios, estimates), display=c("s", "s", "s", rep("e", 7)), caption = "Unreliability estimates for PMC and turnip-style estimation methods. "), include.rownames=F, print.result=F)
estimatesXTable <- gsub("([[:digit:]]+\\.[[:digit:]]*)e(-?)([[:digit:]]+)", "$\\1 \\\\times 10^\\{\\2\\3}$", estimatesXTable)
estimatesXTable <- gsub("0.00e+00", "0", estimatesXTable, fixed=T)

timesXTable <- print(xtable(cbind(scenarios, timesPerSample), display=c("s", "s", "s", rep("e", 7)), caption = "Seconds per sample for different PMC and turnip-style estimation methods. "), include.rownames=F, print.result=F)
timesXTable <- gsub("([[:digit:]]+\\.[[:digit:]]*)e(-?)([[:digit:]]+)", "$\\1 \\\\times 10^\\{\\2\\3}$", timesXTable)
timesXTable <- gsub("0.00e+00", "0", timesXTable, fixed=T)

workNormalizedXTable <- print(xtable(cbind(scenarios, relativeErrors), display=c("s", "s", "s", rep("e", 7)), caption = "Work normalized relative error for PMC and turnip-style estimation methods. "), include.rownames=F, print.result=F)
workNormalizedXTable <- gsub("([[:digit:]]+\\.[[:digit:]]*)e(-?)([[:digit:]]+)", "$\\1 \\\\times 10^\\{\\2\\3}$", workNormalizedXTable)

cat("\\documentclass{article}\\usepackage[landscape]{geometry}\\begin{document}", workNormalizedXTable, estimatesXTable, timesXTable, "\\end{document}", file = "tables.tex")
system("pdflatex tables.tex")
unlink("tables.aux")
unlink("tables.log")
unlink("tables.tex")