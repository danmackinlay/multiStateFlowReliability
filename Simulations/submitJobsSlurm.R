source("./generateScenarios.R")
if(!file.exists("./results")) dir.create("results")
if(!exists("maxJobs")) maxJobs <- nrow(scenarios)
submittedJobs <- 0
for(i in 1:nrow(scenarios))
{
	submit <- TRUE
	resultFile <- file.path("results", scenarios[i, "file"])
	if(!file.exists(resultFile))
	{
		submit <- TRUE
	}
	else
	{
		rm(results)
		load(resultFile)
		if(length(results) != scenarios[i, "nReps"])
		{
			submit <- TRUE
		}
		else submit <- FALSE
		cat(i, ", ", length(results), " / ", scenarios[i, "nReps"],"\n")
	}
	if(submit)
	{
		system2(command = "sbatch", args = c(paste0("--export=SCENARIO_INDEX=", i), "submitScriptSlurm.sh"), wait=TRUE)
		submittedJobs <- submittedJobs + 1
		if(submittedJobs == maxJobs) break
	}
}
