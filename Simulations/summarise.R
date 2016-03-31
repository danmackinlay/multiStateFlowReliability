source("./generateScenarios.R")
allResults <- list()
for(i in 1:nrow(scenarios))
{
	path <- file.path("results", scenarios[i, "file"])
	if(file.exists(path))
	{
		load(path)
		allResults[[i]] <- results
	}
}
secondsPerRun <- lapply(allResults, function(x) unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
averageSecondsPerRun <- unlist(lapply(secondsPerRun, function(x) if(is.null(x)) NA else mean(x)))
secondsSingleSample <- unlist(lapply(allResults, function(x)
	{
		if(is.null(x) || class(x[[1]]) == "generalisedSplittingResult") 
		{
			return(NA)
		}
		else 
		{
			return(mean(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs")/y@n))))
		}
	}))

averageEstimatesFunc <- function(x)
{
	if(is.null(x))
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult") 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data)))))
	}
	else if(class(x[[1]]) == "generalisedSplittingResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimateFirstMoment)))))
	}
}
averageEstimates <- do.call(c, lapply(allResults, averageEstimatesFunc))

varianceSingleSampleFunc <- function(x)
{
	if(is.null(x))
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data*(1 - y@data))))))
	}
	else if(class(x[[1]]) == "generalisedSplittingResult")
	{
		return(NA)
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@varianceEstimate)))))
	}
}
varianceSingleSample <- do.call(c, lapply(allResults, varianceSingleSampleFunc))
workNormalizedSingleSampleVariance <- as.numeric(varianceSingleSample * secondsSingleSample)

varianceFunc <- function(x)
{
	if(is.null(x))
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data*(1 - y@data)/y@n)))))
	}
	else if(class(x[[1]]) == "generalisedSplittingResult")
	{
		return(as.numeric(var(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@varianceEstimate/y@n)))))
	}
}
variances <- do.call(c, lapply(allResults, varianceFunc))
workNormalizedVariance <- as.numeric(variances * averageSecondsPerRun)

save(allResults, secondsSingleSample, varianceSingleSample, workNormalizedSingleSampleVariance, averageEstimates, averageSecondsPerRun, secondsPerRun, variances, workNormalizedVariance, file = "summarised.RData")
