source("./generateScenarios.R")
library(Rmpfr)
allResults <- pilots <- vector(mode = "list", length = nrow(scenarios))
for(i in 1:nrow(scenarios))
{
	path <- file.path("results", scenarios[i, "file"])
	if(file.exists(path))
	{
		load(path)
		allResults[[i]] <- results
		if(scenarios[i, "method"] == "gsFS")
		{
			pilots[[i]] <- pilot
		}
	}
}
averageSecondsPerRun <- unlist(lapply(as.list(1:nrow(scenarios)), function(index)
	{
		x <- allResults[[index]]
		if(is.null(x))
		{
			NA
		}
		else
		{
			mean(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
		}
	}))
averageSecondsPlusPilotPerRun <- unlist(lapply(as.list(1:nrow(scenarios)), function(index)
	{
		x <- allResults[[index]]
		if(is.null(x) || length(x) == 0)
		{
			return(NA)
		}
		else if(class(x[[1]]) == "generalisedSplittingFixedFactorsEvolutionResult")
		{
			pilotTime <- difftime(pilots[[index]]@end, pilots[[index]]@start, units = "secs")
			sampleTime <- mean(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
			return(sampleTime + pilotTime)
		}
		else
		{
			sampleTime <- mean(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
			return(sampleTime)
		}
	}))
averageEstimatesFunc <- function(x)
{
	if(is.null(x) || length(x) == 0)
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult") 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data)))))
	}
	else if(class(x[[1]]) %in% "generalisedSplittingFixedFactorsEvolutionResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@firstMomentSingleSample)))))
	}
}
averageEstimates <- do.call(c, lapply(allResults, averageEstimatesFunc))

#This variance computation just uses the estimates of each of the 100 runs. 
varianceFunc <- function(x)
{
	if(is.null(x) || length(x) == 0)
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data*(1 - y@data)/y@n)))))
	}
	else if(class(x[[1]]) %in% "generalisedSplittingFixedFactorsEvolutionResult")
	{
		values <- do.call(c, lapply(x, function(y) y@estimate))
	}
	else 
	{
		values <- do.call(c, lapply(x, function(y) y@firstMomentSingleSample))
	}
	values <- mpfr(values, prec = 3*getPrec(values))
	sum <- sum(values)
	sumSquared <- sum(values*values)
	return(as.numeric(sumSquared / length(x) - (sum / length(x))^2))
}
variances <- do.call(c, lapply(allResults, varianceFunc))

workNormalizedRelativeVariancePlusPilot <- as.numeric(variances * averageSecondsPlusPilotPerRun / (averageEstimates^2))
relativeErrors <- sqrt(variances) / averageEstimates

save(allResults, averageSecondsPlusPilotPerRun, workNormalizedRelativeVariancePlusPilot, averageEstimates, averageSecondsPerRun, variances, relativeErrors, file = "summarised.RData")
