source("./generateScenarios.R")
if(!file.exists("./results")) dir.create("results")
if(!exists("maxJobs")) maxJobs <- nrow(scenarios)
submittedJobs <- 0
for(i in 1:nrow(scenarios))
{
	resultFile <- file.path("results", scenarios[i, "file"])
	if(!file.exists(resultFile))
	{
		submit <- TRUE
	}
	else
	{
		rm(results)
		load(resultFile)
		if(length(results) != 100)
		{
			submit <- TRUE
		}
		else submit <- FALSE
	}
	if(submit)
	{
		system2(command = "qsub", args = "submitScriptTorque.sh", env = paste0("SCENARIO_INDEX=", i), wait=TRUE)
		submittedJobs <- submittedJobs + 1
		if(submittedJobs == maxJobs) break
	}
}
