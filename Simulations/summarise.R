source("./generateScenarios.R")
results <- list()
for(i in 1:nrow(scenarios))
{
	load(file.path("results", scenarios[i, "file"]))
	results[[i]] <- result
}
times <- unlist(lapply(results, function(x) difftime(x@end, x@start, units = "secs")))
timeSingleSample <- times / unlist(lapply(results, function(x) x@n))

averageEstimatesFunc <- function(x)
{
	if(class(x) == "crudeMCResult") x@data
	else x@estimateFirstMoment
}
averageEstimates <- do.call(c, lapply(results, averageEstimatesFunc))

varianceSingleSampleFunc <- function(x)
{
	if(class(x) == "crudeMCResult") x@data*(1 - x@data)
	else x@varianceEstimate
}
varianceSingleSample <- do.call(c, lapply(results, varianceSingleSampleFunc))

workNormalizedSingleSampleVariance <- as.numeric(varianceSingleSample * timeSingleSample)
averageEstimates <- as.numeric(averageEstimates)
varianceSingleSample <- as.numeric(varianceSingleSample)

save(results, timeSingleSample, varianceSingleSample, workNormalizedSingleSampleVariance, averageEstimates, times, file = "summarised.RData")
