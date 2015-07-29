#Fixed parameters
graphFile <- "../../networkReliability/dodecahedronReordered.graphml"
nCrudeMC <- as.integer(10000000)
nTurnip <- as.integer(1000000)
nTurnipFull <- as.integer(100000)
nPMC <- as.integer(1000000)
nPMCWithTurnip <- as.integer(1000000)
interestVertices <- "0 19"

#Parameters that are varied
distributions <- c(10, 100, 1000, 10000)
thresholds <- 2:4

#crudeMC function
crudeMC <- function(distribution, seed, threshold)
{
	commandLine <- paste0("../x64/Release/crudeMC.exe --threshold ", threshold, " --seed ", seed, " --uniformCapacity 0 ", distribution-1 , " --n ", nCrudeMC, " --graphFile ", graphFile, " --interestVertices ", interestVertices)
	start <- Sys.time()
	result <- system(commandLine, intern=TRUE)
	end <- Sys.time()
	if(substr(result, 1, 39) != "Unreliability probability estimate was ")
	{
		stop("Unexpected output from crudeMC.exe")
	}
	estimate <- as.double(substr(result, 40, nchar(result)))
	relativeError <- sqrt(estimate * (1 - estimate) / nCrudeMC) / estimate
	return(list(time = unclass(difftime(end, start, units="mins"))[1], relativeError = relativeError, estimate = estimate))
}
#PMC function
pmc <- function(distribution, seed, threshold, turnip = FALSE)
{
	commandLine <- paste0("../x64/Release/PMC.exe --threshold ", threshold, " --seed ", seed, " --uniformCapacity 0 ", distribution-1, " --n ", nPMC, " --graphFile ", graphFile, " --interestVertices ", interestVertices)
	if(turnip)
	{
		commandLine <- paste0(commandLine, " --turnip")
	}
	start <- Sys.time()
	result <- system(commandLine, intern=TRUE)
	end <- Sys.time()
	if(substr(result[1], 1, 39) != "Unreliability probability estimate was ")
	{
		stop("Unexpected output from PMC.exe")
	}
	estimate <- as.double(substr(result[1], 40, nchar(result[1])))
	if(substr(result[2], 1, 19) != "Relative error was ")
	{
		stop("Unexpected output from PMC.exe")
	}
	relativeError <- as.numeric(substr(result[2], 20, nchar(result[2])))
	return(list(time = unclass(difftime(end, start, units="mins"))[1], relativeError = relativeError, estimate = estimate))
}

turnip <- function(distribution, seed, threshold, turnip = FALSE, full = FALSE, fullIncrement = 1)
{
	if(full) nSamples <- nTurnipFull
	else nSamples <- nTurnip
	commandLine <- paste0("../x64/Release/turnip.exe --threshold ", threshold, " --seed ", seed, " --uniformCapacity 0 ", distribution-1, " --n ", nTurnip, " --graphFile ", graphFile, " --interestVertices ", interestVertices)
	if(full)
	{
		commandLine <- paste0(commandLine, " --full ", fullIncrement)
	}
	start <- Sys.time()
	result <- system(commandLine, intern=TRUE)
	end <- Sys.time()
	if(substr(result[1], 1, 39) != "Unreliability probability estimate was ")
	{
		stop("Unexpected output from Turnip.exe")
	}
	estimate <- as.double(substr(result[1], 40, nchar(result[1])))
	if(substr(result[2], 1, 19) != "Relative error was ")
	{
		stop("Unexpected output from Turnip.exe")
	}
	relativeError <- as.numeric(substr(result[2], 20, nchar(result[2])))
	return(list(time = unclass(difftime(end, start, units="mins"))[1], relativeError = relativeError, estimate = estimate))
}

seed <- 1
scenarios <- expand.grid(distributions, thresholds)
colnames(scenarios) <- c("Distribution", "Threshold")

#Output data
times <- matrix(0, nrow = nrow(scenarios), ncol = 7)
colnames(times) <- c("crudeMC", "PMC", "PMCPartialTurnip", "TurnipSingleMaxFlow", "TurnipFull1", "TurnipFull2", "TurnipFull3")
estimates <- relativeErrors <- times

#Check directories
if(!file.exists("crudeMC")) dir.create("crudeMC")
if(!file.exists("PMC")) dir.create("PMC")
if(!file.exists("PMCPartialTurnip")) dir.create("PMCPartialTurnip")
if(!file.exists("TurnipSingleMaxFlow")) dir.create("TurnipSingleMaxFlow")
if(!file.exists("TurnipFull1")) dir.create("TurnipFull1")
if(!file.exists("TurnipFull2")) dir.create("TurnipFull2")
if(!file.exists("TurnipFull3")) dir.create("TurnipFull3")

turnipSpecificIncrement <- function(increment)
{
	name <- paste0("TurnipFull", increment)
	turnipFullFile <- file.path(name, paste0(scenarioIndex, ".RData"))
	if(!file.exists(turnipFullFile))
	{
			turnipFullResult <- turnip(distribution = distribution, seed = seed, threshold = threshold, full = TRUE, fullIncrement = increment)
			save(turnipFullResult, file = turnipFullFile)
	}
	else load(turnipFullFile)
	times[scenarioIndex, name] <<- turnipFullResult$time
	relativeErrors[scenarioIndex, name] <<- turnipFullResult$relativeError
	estimates[scenarioIndex, name] <<- turnipFullResult$estimate
	seed <<- seed + 1
}
for(scenarioIndex in 1:nrow(scenarios))
{
	distribution <- scenarios[scenarioIndex, "Distribution"]
	threshold <- scenarios[scenarioIndex, "Threshold"]

	#crudeMC
	crudeMCFile <- file.path("crudeMC", paste0(scenarioIndex, ".RData"))
	if(!file.exists(crudeMCFile))
	{
		crudeMCResult <- crudeMC(distribution = distribution, seed = seed, threshold = threshold)
		save(crudeMCResult, file = crudeMCFile)
	}
	else
	{
		load(crudeMCFile)
	}
	times[scenarioIndex, "crudeMC"] <- crudeMCResult$time
	relativeErrors[scenarioIndex, "crudeMC"] <- crudeMCResult$relativeError
	estimates[scenarioIndex, "crudeMC"] <- crudeMCResult$estimate
	seed <- seed + 1

	#PMC
	pmcFile <- file.path("PMC", paste0(scenarioIndex, ".RData"))
	if(!file.exists(pmcFile))
	{
			pmcResult <- pmc(distribution = distribution, seed = seed, threshold = threshold, turnip=FALSE)
			save(pmcResult, file = pmcFile)
	}
	else load(pmcFile)
	times[scenarioIndex, "PMC"] <- pmcResult$time
	relativeErrors[scenarioIndex, "PMC"] <- pmcResult$relativeError
	estimates[scenarioIndex, "PMC"] <- pmcResult$estimate
	seed <- seed + 1

	#PMC with Turnip
	pmcWithTurnipFile <- file.path("PMCPartialTurnip", paste0(scenarioIndex, ".RData"))
	if(!file.exists(pmcWithTurnipFile))
	{
			pmcWithTurnipResult <- pmc(distribution = distribution, seed = seed, threshold = threshold, turnip=TRUE)
			save(pmcWithTurnipResult, file = pmcWithTurnipFile)
	}
	else load(pmcWithTurnipFile)
	times[scenarioIndex, "PMCPartialTurnip"] <- pmcWithTurnipResult$time
	relativeErrors[scenarioIndex, "PMCPartialTurnip"] <- pmcWithTurnipResult$relativeError
	estimates[scenarioIndex, "PMCPartialTurnip"] <- pmcWithTurnipResult$estimate
	seed <- seed + 1

	#Turnip with single edge max flow
	turnipFile <- file.path("TurnipSingleMaxFlow", paste0(scenarioIndex, ".RData"))
	if(!file.exists(turnipFile))
	{
			turnipResult <- turnip(distribution = distribution, seed = seed, threshold = threshold)
			save(turnipResult, file = turnipFile)
	}
	else load(turnipFile)
	times[scenarioIndex, "TurnipSingleMaxFlow"] <- turnipResult$time
	relativeErrors[scenarioIndex, "TurnipSingleMaxFlow"] <- turnipResult$relativeError
	estimates[scenarioIndex, "TurnipSingleMaxFlow"] <- turnipResult$estimate
	seed <- seed + 1

	#Turnip with full
	for(increment in 1:3)
	{
		turnipSpecificIncrement(increment)
	}
}