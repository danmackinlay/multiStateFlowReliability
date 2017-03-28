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
secondsSingleSample <- unlist(lapply(as.list(1:nrow(scenarios)), function(index)
	{
		x <- allResults[[index]]
		if(is.null(x) || length(x) == 0)
		{
			return(NA)
		}
		else 
		{
			return(mean(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs")/y@n))))
		}
	}))
secondsSingleSamplePlusPilot <- unlist(lapply(as.list(1:nrow(scenarios)), function(index)
	{
		x <- allResults[[index]]
		if(is.null(x) || length(x) == 0)
		{
			return(NA)
		}
		else if(class(x[[1]]) == "generalisedSplittingFixedFactorsEvolutionResult")
		{
			nSamples <- length(x) * x[[1]]@n
			pilotTime <- difftime(pilots[[index]]@end, pilots[[index]]@start, units = "secs")
			sampleTime <- sum(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
			return((sampleTime + pilotTime) / nSamples)
		}
		else
		{
			nSamples <- length(x) * x[[1]]@n
			sampleTime <- sum(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
			return(sampleTime/nSamples)
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

varianceSingleSampleFunc <- function(index)
{
	x <- allResults[[index]]
	if(is.null(x) || length(x) == 0)
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data*(1 - y@data))))))
	}
	else if(class(x[[1]]) %in% "generalisedSplittingFixedFactorsEvolutionResult")
	{
		return(variances[index]*x[[1]]@n)
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@varianceSingleSample)))))
	}
}
varianceSingleSample <- do.call(c, lapply(1:length(allResults), varianceSingleSampleFunc))
workNormalizedSingleSampleVariancePlusPilot <- as.numeric(varianceSingleSample * secondsSingleSamplePlusPilot)
workNormalizedSingleSampleVariance <- as.numeric(varianceSingleSample * secondsSingleSample)
relativeErrors <- sqrt(variances) / averageEstimates

save(allResults, secondsSingleSample, varianceSingleSample, workNormalizedSingleSampleVariancePlusPilot, workNormalizedSingleSampleVariance, averageEstimates, averageSecondsPerRun, variances, relativeErrors, secondsSingleSamplePlusPilot, file = "summarised.RData")
