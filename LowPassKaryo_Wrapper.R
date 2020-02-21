
##This will need to be customised to the site:
SITE.CONFIG.FILE <- ##</absolute/path/to/config.R>
source(SITE.CONFIG.FILE)

				###
				###
				###	This front end script is intended to be the only one invoked directly by the end user.
				###
				###


these.args <- commandArgs(trailing=TRUE)
	if( DEBUG )
	{
		if( length(these.args) == 0 )
		{
		these.args <- c("--design", "test_design_file.txt")
		}
	}


		########
		##MAIN##
		########
########################################################################################################################
main <- function(
		ARGS
		)
{
init.dir <- getwd()
on.exit(setwd(init.dir))


	##
	##Parse the command line arguments:
	###################################

design.file.name <- parse.command.line(passed.args = ARGS)

	if( is.na(design.file.name))
	{
	return()
	}

GENOME.LOOKUP.TABLE <- as.matrix(read.delim(GENOMES.LOOKUP.FILE, sep="\t", header=TRUE, comment = "", quote=""))
PARAMS <- parse.design.file( DESIGN.FILE.NAME = design.file.name, GENOME.LOOKUP= GENOME.LOOKUP.TABLE )

	##
	##Next up: associate full fastq paths with the LIMS IDs and determine whether it's a single or paired end run.
	##############################################################################################################


FASTQ.FILES <- associate.fastq.with.lims.id(
					params=PARAMS
					)
row.names(FASTQ.FILES) <- PARAMS[["DESIGN"]][, "LIMS_ID"]

PARAMS[["FASTQ_FILES"]] <- FASTQ.FILES


these.genome.tags <- PARAMS[["DESIGN"]][, "GENOME_TAG"]

GENOME.FILES <- rep(NA, length(these.genome.tags))
QDNA.FILES <- rep(NA, length(these.genome.tags))

	for( i in seq_along(these.genome.tags) )
	{
	this.tag <- these.genome.tags[i]

	w <- which( GENOME.LOOKUP.TABLE[ ,"Tag"] == this.tag )
		if( length(w) != 1 )
		{
		msg <- paste("Out of Coffee Error: could not find a unique match to the tag \"", this.tag, "\" in the genomes lookup table.\n", 
					"This should be impossible at this stage.\n", sep="")
		stop(msg)
		}
	this.reference <- GENOME.LOOKUP.TABLE[ , "Reference_Path"][w]
	this.qdna.path <- GENOME.LOOKUP.TABLE[, "QDNAseq_Annotated_Bins_Path"][w]

	GENOME.FILES[i] <- this.reference
	QDNA.FILES[i] <- this.qdna.path
	}


PARAMS[["GENOME_FILES"]] <- GENOME.FILES
PARAMS[["QDNA_FILES"]] <- QDNA.FILES


DETAILS <- cbind( FASTQ.FILES, GENOME.FILES, QDNA.FILES)
colnames(DETAILS) <- c( "R1", "R2", "Ref", "Qdna")

	###########################################################
	###							###
	###	Generate Basic Output Directory Structure	###
	###							###
	###########################################################


##
##set up project dir
##raw fastq dir and soft links
##FASTQC_RAW
##FASTQ_TRIMMED
##FASTQC_TRIMMED
##TMP dir
##SAM
##BAM
##BAM_SORTED
##REPORTS
##



	if( !file.exists(PARAMS[[ "OUTPUT_ROOT_DIR" ]] ))
	{
	cat("Attempting to create output root directory:\n")

	cmd <- paste("mkdir ", PARAMS[["OUTPUT_ROOT_DIR"]], sep="")
	S <- system(cmd, intern = FALSE)
		if( S != 0)
		{
		msg <- paste("Failed to create output root directory: please check the path:\n\t\"", PARAMS[["OUTPUT.ROOT.DIRECTORY"]], "\"\n", sep="")
		stop(msg)
		}

	}

cat("Checking project name for suitability as a directory name:\n")
project.base <- sanitise.names(PARAMS[["PROJECT_NAME"]])$clean
PARAMS[["PROJECT_ROOT"]] <- project.base

PARAMS[["WORKING_DIR"]] <- paste(PARAMS[["OUTPUT_ROOT_DIR"]], PARAMS[["PROJECT_ROOT"]], "/", sep="")

	if( !file.exists(PARAMS[["WORKING_DIR"]]) )
	{
	cat("Attempting to create output directory for this project.\n")
	cmd <-paste("mkdir ", PARAMS[["WORKING_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)

		if( S != 0 )
		{
		msg <- paste("Failed to create project working directory:\n\t\"", PARAMS[["WORKING_DIR"]], "\"\n", sep="")
		stop(msg)
		}
	}


	###########################################################################
	###	Create FASTQ_RAW directory and soft link the raw fastq files	###
	###########################################################################

setwd(PARAMS[["WORKING_DIR"]])

PARAMS[["FASTQ_RAW_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "FASTQ_RAW", sep=""))
	if( !file.exists(PARAMS[["FASTQ_RAW_DIR"]]) )
	{
	cmd <- paste("mkdir ", PARAMS[["FASTQ_RAW_DIR"]], sep="")
	cat("Attempting to create local FastQ directory.\n")

	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Failed to create local fastq directory:\n\t\"", PARAMS[["FASTQ_RAW_DIR"]], "\"\n", sep="")
		stop(msg)
		}
	}



setwd(PARAMS[["FASTQ_RAW_DIR"]])

LIMS.IDS <- PARAMS[["DESIGN"]][, "LIMS_ID"]
SAMPLE.NAMES.RAW <- PARAMS[["DESIGN"]][, "SAMPLE_NAME_RAW"]

	for(i in seq_along(LIMS.IDS) )
	{
	this.lims.id <- LIMS.IDS[i]
	these.fastq.files <- PARAMS[["FASTQ_FILES"]][ this.lims.id,,drop=FALSE]

	
		for(j in seq_along(these.fastq.files))
		{
		this.source.fastq.file <- these.fastq.files[j]
			if( this.source.fastq.file == "SINGLE_END" )
			{

			next()
			}
		full.source.fastq.file <- paste(PARAMS[[ "FASTQ_DIR" ]], this.source.fastq.file, sep="")

		this.local.fastq.file <- this.source.fastq.file

			if( !file.exists( this.local.fastq.file ) )
			{

				if( !file.exists(full.source.fastq.file) )
				{
				msg <- paste("Error locating orig fastq files: could not find the following expected file:\n\t\"", full.source.fastq.file, "\"\n", sep="")
				stop(msg)
				}

			ln.cmd <- paste("ln -s ", full.source.fastq.file, " ", this.local.fastq.file, sep="")
			S <- system(ln.cmd, intern = FALSE )
				if( S != 0 )
				{
				msg <- paste("Softlinking failed for \"", this.source.fastq.file, "\"\n", sep="")
				stop(msg)
				}
			}
		}
	}


	##
	##Create the individual output directories
	##

PARAMS[["FASTQC_RAW_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "FASTQC_RAW", sep=""))
	if( !file.exists( PARAMS[["FASTQC_RAW_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["FASTQC_RAW_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(length(S) > 0 & S != 0 )##not entirely sure why this particular one is not returning anything on success
		{
		stop("Attempt to create FastQC_Raw output directory failed.\n")
		} 
	}


PARAMS[["FASTQ_TRIMMED_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "FASTQ_TRIMMED", sep=""))
	if( !file.exists( PARAMS[["FASTQ_TRIMMED_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["FASTQ_TRIMMED_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create FastQ_Trimmed output directory failed.\n")
		} 
	}


PARAMS[["FASTQC_TRIMMED_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "FASTQC_TRIMMED", sep=""))
	if( !file.exists( PARAMS[["FASTQC_TRIMMED_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["FASTQC_TRIMMED_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create FastQC_Trimmed output directory failed.\n")
		} 
	} 


PARAMS[["TMP_DIR"]] <- path.cap("TMP")	##TODO: add an option to sit this somewhere other than within output dir - like on scratch
	if(nchar(PARAMS[["TMP_DIR"]])>0)
	{
	full.tmp <- paste(PARAMS[["WORKING_DIR"]], PARAMS[["TMP_DIR"]], sep="")
		if( !file.exists( full.tmp ) )
		{
		cmd <- paste("mkdir ", full.tmp, sep="")
		S <- system(cmd, intern=FALSE)
			if(S != 0 )
			{
			stop("Attempt to create TMP_DIR output directory failed.\n")
			} 
		} 	

	}	

PARAMS[["SAM_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], PARAMS[["TMP_DIR"]], "SAM", sep=""))
	if( !file.exists( PARAMS[["SAM_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["SAM_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create SAM_DIR output directory failed.\n")
		} 
	} 




PARAMS[["BAM_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], PARAMS[["TMP_DIR"]], "BAM", sep=""))
	if( !file.exists( PARAMS[["BAM_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["BAM_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create BAM_DIR output directory failed.\n")
		} 
	} 


PARAMS[["BAM_SORTED_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "BAM_SORTED", sep=""))
	if( !file.exists( PARAMS[["BAM_SORTED_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["BAM_SORTED_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create BAM_SORTED_DIR output directory failed.\n")
		} 
	} 


PARAMS[["QDNASEQ_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "QDNASEQ", sep=""))
	if( !file.exists( PARAMS[["QDNASEQ_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["QDNASEQ_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create QDNASEQ_DIR output directory failed.\n")
		} 
	} 

PARAMS[["XY_QDNASEQ_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "QDNASEQ_XY", sep=""))
	if( !file.exists( PARAMS[["XY_QDNASEQ_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["XY_QDNASEQ_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create XY_QDNASEQ_DIR output directory failed.\n")
		} 
	} 




PARAMS[["RESULTS_DIR"]] <- path.cap(paste(PARAMS[["WORKING_DIR"]], "RESULTS", sep=""))
	if( !file.exists( PARAMS[["RESULTS_DIR"]] ) )
	{
	cmd <- paste("mkdir ", PARAMS[["RESULTS_DIR"]], sep="")
	S <- system(cmd, intern=FALSE)
		if(S != 0 )
		{
		stop("Attempt to create RESULTS_DIR output directory failed.\n")
		} 
	} 



##save the translations as an R object so that the summary builder can use them later.

params.lookup.file <- paste(PARAMS[["RESULTS_DIR"]], PARAMS[["PROJECT_NAME"]], "_PARAMS.RDat", sep="")
PARAMS[["PARAMS_LOOKUP_FILE"]] <- params.lookup.file
save(PARAMS, file=params.lookup.file)




	###################################################
	###						###
	###	Set Up the Scheduler Submissions	###
	###						###
	###################################################

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


DEPENDENCIES <- vector()



FORCE.SUMMARY.RUN <- FALSE	##if this is set to TRUE, then the summary job will be submitted even if the completed tag file exists.
				##this wil then enable an automatic re-run of the summary is we update or add any individual samples


loc.max.threads <- DEFAULT.MAX.THREADS
##check to see if we can find a cpu setting in PARAMS[["SCHEDULER_OPTIONS"]]

tmp <- unlist(strsplit(gsub(" +", " ", PARAMS[["SCHEDULER_OPTIONS"]]), split=" "))

g <- grep(THREADS.REGEXP, tmp)
	if(length(g) == 1 )
	{

	TMP <- as.numeric(gsub(" +", "", gsub(THREADS.REGEXP, "", tmp[g])))
	
		if( is.na(TMP) )
		{
		msg <- paste("Out of Coffee Error: could not infer a sensible number of CPUs from the following scheduler options:\n\t\"",
					PARAMS[[ "SCHEDULER_OPTIONS" ]], "\"\n", sep="")

		msg <- paste(msg, "\ntmp[", g, "]==\"", tmp[g], "\"\n", sep="")
		stop(msg)			
		}

	loc.max.threads <- TMP
	}

	if( SUBMIT.FROM.ENGINE )
	{##give up a core for running the instance

	loc.max.threads <- loc.max.threads -1
	}


	for(i in seq_along(LIMS.IDS))
	{
	this.lims.id <- LIMS.IDS[i]
	this.sample.name.raw <- SAMPLE.NAMES.RAW[i]
	this.completed.file <- paste(PARAMS[["RESULTS_DIR"]], ".", this.lims.id, "_individual.done", sep="")
	
		if( TRUE || !file.exists(this.completed.file) )
		{

		this.log.file <- paste(PARAMS[["WORKING_DIR"]], "Engine_", this.lims.id, "_log.txt", sep="")
		re.run.ind <- 0
			while( file.exists(this.log.file) )
			{
			re.run.ind <- re.run.ind + 1

			this.log.file <- paste(PARAMS[["WORKING_DIR"]], "Engine_", this.lims.id, "_restart-", re.run.ind, "_log.txt", sep="")
			}

		FQ2 <- PARAMS[["FASTQ_FILES"]][i,2]
			if( FQ2 != "SINGLE_END" )
			{
			FQ2 <- paste(PARAMS[["FASTQ_RAW_DIR"]], FQ2, sep="")

			}

		these.params <- paste(
				"\"CONFIG.R_PATH=", SITE.CONFIG.FILE, "\" ",
				"\"SCHEDULER_OPTIONS=", PARAMS[["SCHEDULER_OPTIONS"]], "\" ",
				"\"MAX_THREADS=", loc.max.threads, "\" ",	##this is really only used at the alignment stage, everything else runs single threaded.
				"\"PROJECT_NAME=", PARAMS[["PROJECT_NAME"]], "\" ",
				"\"PROJECT_DIR=",PARAMS[["WORKING_DIR"]], "\" ",
				"\"LIMS_ID=", this.lims.id, "\" ",
				"\"SAMPLE_NAME_RAW=", this.sample.name.raw, "\" ",
				"\"FASTQ_1=", PARAMS[["FASTQ_RAW_DIR"]], PARAMS[["FASTQ_FILES"]][i,1], "\" ",
				"\"FASTQ_2=", FQ2, "\" ",
				"\"REFERENCE_PATH=", PARAMS[["GENOME_FILES"]][i], "\" ",
				"\"ANNOTATED_BINS_FILE=",	PARAMS[["QDNA_FILES"]][i], "\" ",
				"\"FASTQC_RAW_DIR=",PARAMS[["FASTQC_RAW_DIR"]], "\" ",
				"\"FASTQ_TRIMMED_DIR=", PARAMS[["FASTQ_TRIMMED_DIR"]], "\" ",
				"\"FASTQC_TRIMMED_DIR=", PARAMS[["FASTQC_TRIMMED_DIR"]], "\" ",
				"\"SAM_DIR=", PARAMS[["SAM_DIR"]], "\" ",
				"\"BAM_DIR=", PARAMS[["BAM_DIR"]], "\" ",
				"\"BAM_SORTED_DIR=",PARAMS[["BAM_SORTED_DIR"]], "\" ",
				"\"QDNASEQ_BED_FILE=", PARAMS[["DESIGN"]][,"ANNO_BED"][i], "\" ",
				"\"QDNASEQ_DIR=", PARAMS[["QDNASEQ_DIR"]], "\" ",
				"\"XY_QDNASEQ_DIR=",  PARAMS[["XY_QDNASEQ_DIR"]], "\" ",
				"\"RESULTS_DIR=", PARAMS[["RESULTS_DIR"]], "\" ",
				"\"COMPLETION_FILE_NAME=", this.completed.file, "\" ",
				sep="")

			###Interactive session version
			#CMD <- paste(RSCRIPT, " ", SV_ENGINE.R, " ", these.params, sep="")
			#S <- system(CMD, intern=FALSE)
			###############################


		
		this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "SV_ENGINE.R"]], " ", these.params, sep="")
		this.job.name <- paste(this.lims.id, "_Eng", sep="")


		sub.cmd <- build.scheduler.submission(
							scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
							job.name = this.job.name,
							log.file = this.log.file,
							job.cmd = this.job.cmd
							)



		S <- system(sub.cmd, intern=TRUE)
		g <- grep("^Submitted batch job ", S)
			if( length(g) == 0 )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n", sep="")
			stop(msg)
			}

		msg <- paste("Job \"", this.lims.id, "\" submitted:\n\t\"", S, "\"\n", sep="")
		cat(msg)
		FORCE.SUMMARY.RUN <- TRUE
		this.id <- get.scheduler.id(S)  


		DEPENDENCIES <- c( DEPENDENCIES, this.id)

		}
	}##end invididual samples processing loop.


	##
	##Set up and submit the summarisation job which will merge the results from all of the individual samples
	##

completion.file.name <- paste(PARAMS[["RESULTS_DIR"]],COMPLETION.MARKER.FILE.NAME , sep="")

	if( FORCE.SUMMARY.RUN || !file.exists( completion.file.name ) )
	{

		if( FORCE.SUMMARY.RUN && file.exists(completion.file.name) )
		{
		cat("Re-running summary construction due to re-run of one or more individual samples.\n")
		}

	deps <- NULL
		if( length(DEPENDENCIES) > 0 )
		{
		deps <- DEPENDENCIES
		}
	these.params <- paste(
				"\"CONFIG.R_PATH=", SITE.CONFIG.FILE, "\" ",	
				"\"PARAMS_LOOKUP_FILE=", PARAMS[["PARAMS_LOOKUP_FILE"]], "\" ",
				"\"COMPLETION_FILE_NAME=", completion.file.name, "\" ",
				sep="")

	this.job.cmd <- paste(
				BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "PREP_REPORT.R" ]] , " ",
						these.params, 
						sep="")

	this.job.name <- paste(PARAMS[["PROJECT_NAME"]], "_Summary", sep="")
	this.log.file <- paste(PARAMS[["WORKING_DIR"]], "Project_", PARAMS[["PROJECT_NAME"]], "_Summary_Job_log.txt", sep="")

	sub.cmd <-  build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)


	msg <- paste("\n\nSubmitting summarisation job: on (successful) completion, this will write:\n\t\"", completion.file.name, "\"\n\n", sep="")
	cat(msg)

	S <- system(sub.cmd, intern=TRUE)
	this.id <- get.scheduler.id(S) 
		if( is.na(this.id) )
		{
		msg <- paste("Error submitting summary job\n\t\"", S, "\"\n", sep="")
		stop(msg)
		}
 	msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
	cat(msg)

	}



return(0)
}##end MAIN
########################################################################################################################


a <- main(ARGS = these.args)
