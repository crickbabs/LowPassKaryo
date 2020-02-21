


	##################################################################################################
	##												##
	##	If you use a scheduler other than slurm, please add replacement submisison functions	##
	##		here and update the variable settings at the bottom of the file			##	
	##												##
	##################################################################################################



	###################################
	###				###
###########	SCHEDULER: SLURM	###########
	###				###
	###################################


	##########################
	##BUILD.SLURM.SUBMISSION##
	##########################
########################################################################################################################
build.slurm.submission <- function(
			scheduler.options = " --ntasks=1 -t 3-00:00:00 ",
			dependencies = NULL,
			job.name,
			log.file,
			job.cmd
			)
{



cat("Constructing slurm submission statement.\n")
sub.cmd <- paste("sbatch ",
				" ", scheduler.options, " ",
				" --job-name ", job.name, " ",
				" -o ", log.file,
				sep="")
	if( !is.null(dependencies) )
	{

	sub.cmd <- paste(sub.cmd, 
		" --depend=afterok:", paste(dependencies, collapse=":"), " ", 
		" --kill-on-invalid-dep=TRUE " , sep="")
	}
sub.cmd <- paste(sub.cmd, 
				" --wrap='", job.cmd, "'",
				sep="")



return(sub.cmd)
}##end BUILD.SLURM.SUBMISSION
########################################################################################################################


	################
	##GET.SLURM.ID##
	################
#
#Extract the scheduler Job ID from the return message
#
########################################################################################################################
get.slurm.id <- function(
			X
			)
{

ret <- unlist(strsplit(X, split=" "))[4]

return(ret)
}##end GET.SLURM.ID
########################################################################################################################



	######################
	##CHECK.SLURM.FOR.ID##
	######################
#
#how do we check the queue for a particular job id?
#
########################################################################################################################
check.slurm.for.id <- function(
			ID
			)
{


cmd <- paste("squeue | grep ", ID, sep="")
S <- suppressWarnings(system(cmd, intern=TRUE))	##if the id isn't found, system will throw a warning that grep failed. 
						##	Since we're explicitly dealign with the case were no result is found
						##	below the warning is both unnecessary and potentiall confusing in the 
						##	log file, so we squelch it using supressWarnings()
						##
	if( length(S) > 0 )
	{
	ret <- S
	}else
	{
	ret <- NA
	}
return(ret)

}##end CHECK.SLURM.FOR.ID
########################################################################################################################




	###########################################
	###					###
###########	SET SCHEDULER FUNCTIONS		###########
	###					###
	###########################################

build.scheduler.submission <- build.slurm.submission	##if you don't use slurm as a job scheduler then update these assignments!
get.scheduler.id <- get.slurm.id
check.queue.for.id <- check.slurm.for.id
DEFAULT.SCHEDULER.OPTIONS <- paste(
					" --ntasks=1 ",
					" -t 1-00:00:00 ",	##set a walltime of 1 day.
					" --cpus-per-task=8 ",
					sep="")

THREADS.REGEXP <- "--cpus-per-task="	##a regular expression to match the cpu request in a submission
