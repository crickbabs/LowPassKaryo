

init.dir <- getwd()

	###################################################################
	##this will sit on a node and invoke the various jobs per-patient##
	###################################################################


##
##Dispatch a wake up message
##
host <- system("hostname", intern=TRUE)
msg <- paste("SV_ENGINE.R starting run on \"", host, "\" at ", date(), "\n\n", sep="")
cat(msg)
flush.console()

my.args <- commandArgs(trailingOnly=TRUE)


EXPECTED.PARAMS <- c(		
		"CONFIG.R_PATH",
		"SCHEDULER_OPTIONS",
		"MAX_THREADS",
		"PROJECT_NAME",
		"PROJECT_DIR",
		"LIMS_ID",
		"SAMPLE_NAME_RAW",
		"FASTQ_1",
		"FASTQ_2",
		"REFERENCE_PATH",
		"ANNOTATED_BINS_FILE",	
		"FASTQC_RAW_DIR",
		"FASTQ_TRIMMED_DIR",
		"FASTQC_TRIMMED_DIR",
		"SAM_DIR",
		"BAM_DIR",
		"BAM_SORTED_DIR",
		"QDNASEQ_BED_FILE",
		"QDNASEQ_DIR",
		"XY_QDNASEQ_DIR",
		"RESULTS_DIR",
		"COMPLETION_FILE_NAME"
		)
if( FALSE )
{
##cat("SETTING TESTING PARAMETERS\n")
##my.args <- c()
}


cat("Called with the following parameters:\n\t",
			paste(my.args, collapse="\n\t"),
			"\n\n\n", sep="")

cat(">>>>>", date(), "\tParsing command line params.\n", sep="")


	if( length(my.args) != length(EXPECTED.PARAMS) )
	{
	msg <- paste("This script is intended to be called by another wrapper script, and requires ", length(EXPECTED.PARAMS), " command line parameters (found ", length(my.args), "):\n\t",
			 sep="")

	stop(msg)
	}


candidate <- unlist(strsplit(my.args, split="="))[1]
	if( candidate != EXPECTED.PARAMS[1] )
	{
	msg <- paste("This script is intended to be called by another wrapper script, and must be passed \"CONFIG.R_PATH=\" before any other params.\n")

	stop(msg)	
	}

config.path <- unlist(strsplit(my.args[1], split="="))[2]
	if( !file.exists(config.path) )
	{
	msg <- paste("Config file not found: couldn't read\n\t\"", config.path, "\"\n", sep="")
	stop(msg)
	}

source(config.path)##this gives us, amongst other things, access to vecsplit, which is nice.


tmp <- singlesplit(my.args, split="=")

M <- match( EXPECTED.PARAMS, tmp[,1] )
m <- which(is.na(M))
	if( length(m) > 0 )
	{
	msg <- paste("The following command line parameter names were not detected:\n\t",
				paste(EXPECTED.PARAMS[m], collapse="\n\t"), "\n", sep="")
	stop(msg)
	}


PARAMS <- list()

	for(i in seq_along(EXPECTED.PARAMS))
	{
	this.param <- EXPECTED.PARAMS[i]
	this.value <- tmp[M[i],2]


	PARAMS[[ this.param ]] <- this.value
	}






working.dir <- PARAMS[[ "PROJECT_DIR" ]]
setwd(working.dir)

this.lims.id <- PARAMS[["LIMS_ID"]]
this.sample.name.raw <- PARAMS[["SAMPLE_NAME_RAW"]]

FORCE.PROCESS.DOWNSTREAM <- FALSE 	##if any stage is run then all subsequent steps must be re-run


	###########################
	###			###
	###	FASTQC_RAW	###########################################################################################
	###			###
	###########################
cat(">>>>>", date(), "\tSetting up FastQC call on raw fastq files.\n", sep="")
flush.console()

#################
#EXPECTED.PARAMS <- c(
#		"FASTQC",		##system specific method for invoking FASTQC 
#					##(e.g. including module load statements on camp)
#		"INPUT.FASTQ.FILE",	## fastq file you wish to QC
#		"OUTPUT.DIRECTORY",	##where should FastQC write its results?
#		"COMPLETION.FILE.NAME"	##file to be created to indicate _successful_ completion of the job.
#		)
#################


fastqc.1.raw.done <- paste(PARAMS[["FASTQC_RAW_DIR"]], ".", this.lims.id, "_R1.fastqc_raw.done", sep="")
fastqc.2.raw.done <- paste(PARAMS[["FASTQC_RAW_DIR"]], ".", this.lims.id, "_R2.fastqc_raw.done", sep="")



DEPENDENCIES <- vector()
SUBMISSION.IDS <- vector()



	if( !file.exists(fastqc.1.raw.done) )
	{

	this.input.fastq.file <- PARAMS[["FASTQ_1"]]
	this.output.directory <- PARAMS[["FASTQC_RAW_DIR"]]
	
	these.script.params <- paste(
					"\"FASTQC=", BINARY.CALLS[["FASTQC"]], "\" ",##defined in the config file
					"\"INPUT.FASTQ.FILE=", this.input.fastq.file, "\" ", 
					"\"OUTPUT.DIRECTORY=", this.output.directory, "\" ",
					"\"COMPLETION.FILE.NAME=", fastqc.1.raw.done, "\" ",
					sep="")


	this.job.name <- paste( "FQCr1_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_r1_fqc_log.txt", sep="")

	
	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "FASTQC.R" ]], " ", these.script.params, sep="")



		if( SUBMIT.FROM.ENGINE )
		{
		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
	
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)


		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (raw) for Read 1 with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(sub.cmd, intern=TRUE)
	
		this.id <- get.scheduler.id(S) 
			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}
		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)


		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)

		msg <- paste("Job \"", this.id, "\" submitted:\n\t\"", S, "\"\n", sep="")
		cat(msg)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (raw) for Read 1 directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}


	}else
	{
	msg <- paste("\tFASTQC (raw) for Read 1 already done: skipping!\n")
	cat(msg)
	}


	if( PARAMS[["FASTQ_2"]] != "SINGLE_END" &&  !file.exists(fastqc.2.raw.done) )
	{


	this.input.fastq.file <- PARAMS[["FASTQ_2"]]
	this.output.directory <- PARAMS[["FASTQC_RAW_DIR"]]
	
	these.script.params <- paste(
					"\"FASTQC=", BINARY.CALLS[["FASTQC"]], "\" ",##defined in the config file
					"\"INPUT.FASTQ.FILE=", this.input.fastq.file, "\" ", 
					"\"OUTPUT.DIRECTORY=", this.output.directory, "\" ",
					"\"COMPLETION.FILE.NAME=", fastqc.2.raw.done, "\" ",
					sep="")


	this.job.name <- paste( "FQCr2_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_r2_fqc_log.txt", sep="")
	
	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "FASTQC.R" ]], " ", these.script.params, sep="")

		if( SUBMIT.FROM.ENGINE )
		{
		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
	
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)


		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (raw) for Read 2 with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(sub.cmd, intern=TRUE)
	
		this.id <- get.scheduler.id(S)  
			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}
		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)

		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)

		msg <- paste("Job \"", this.id, "\" submitted:\n\t\"", S, "\"\n", sep="")
		cat(msg)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (raw) for Read 2 directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}


	}else
	{
		if( PARAMS[["FASTQ_2"]] != "SINGLE_END" )
		{
		msg <- paste("\tFASTQC (raw) for Read 2 already done: skipping!\n")
		
		}else
		{
		msg <- paste("\tThis run is single ended so there is no R2 to consider.\n")
		
		}
	cat(msg)
	}




	###########################
	###			###
	###	QCTRIM		###########################################################################################
	###			###
	###########################
cat(">>>>>", date(), "\tSetting up QCTrim call on raw fastq files.\n", sep="")
flush.console()


#################
#EXPECTED.PARAMS <- c(
#               "QCTRIM",           ##system specific method for invoking TRIMGALORE (e.g. including module load statements on camp)#
#		"INPUT.FASTQ.FILE_R1",	##Use this option to pass single end data
#		"INPUT.FASTQ.FILE_R2",	## _IF_ data is single ended, pass NA here _AND_ set the PAIRED.END param to FALSE
#		"PAIRED.END",		##Set TRUE if data is paired end and FALSE if not
#		"OUTPUT.DIRECTORY",
#		"COMPLETION.FILE.NAME"
#		)
#################



trim.raw.done <- paste(PARAMS[["FASTQ_TRIMMED_DIR"]], ".", this.lims.id, ".fastq_trimmed.done", sep="")##R1 & r2 are trimmed as a pair


	if( !file.exists( trim.raw.done ) )##rerunning the QC doesn't automatically mean we need to re-trim so no over-rule here.
	{
	FORCE.PROCESS.DOWNSTREAM <- TRUE

	this.output.directory <- PARAMS[["FASTQ_TRIMMED_DIR"]]

	these.params <- paste(
					"\"QCTRIM=", BINARY.CALLS[["QCTRIM"]], "\" ",
					"INPUT.FASTQ.FILE_R1=", PARAMS[["FASTQ_1"]], " ",
					"INPUT.FASTQ.FILE_R2=", PARAMS[["FASTQ_2"]], " ",
					"PAIRED.END=", (PARAMS[["FASTQ_2"]]!= "SINGLE_END"), " ", 
					"OUTPUT.DIRECTORY=", this.output.directory , " ",
					"COMPLETION.FILE.NAME=", trim.raw.done, " ",
				sep="")

	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[["QCTRIM.R"]] , " ", 
 					these.params,
					sep="")

	

	this.job.name <- paste( "Trim_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_Trim_log.txt", sep="")
	


		if( SUBMIT.FROM.ENGINE )
		{

		deps <- NULL
			if( length(DEPENDENCIES) > 0 )
			{
			deps <- DEPENDENCIES
			DEPENDENCIES <- vector()	##re-set since this is now the only bottle neck
			}

		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)


		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling Trimmer with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(sub.cmd, intern=TRUE)

		this.id <- get.scheduler.id(S) 
 			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}
		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)

		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling QCTrimmer directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}



	}else
	{
	msg <- paste("\tTrimming has already been completed: skipping!\n")
	cat(msg)
	}




	##	
	##At this point we need to recover the names of the _trimmed_ fastq files for this sample
	##

	##
	pre.r1.base <- unlist(strsplit(PARAMS[["FASTQ_1"]], split=.Platform$file.sep))
	pre.r1.base <- pre.r1.base[ length(pre.r1.base) ]

	pre.r2.base <- unlist(strsplit(PARAMS[["FASTQ_2"]], split=.Platform$file.sep))
	pre.r2.base <- pre.r2.base[ length(pre.r2.base) ]

		if(PARAMS[["FASTQ_2"]] != "SINGLE_END" )
		{
		trimmed.r1 <- gsub("\\.fastq\\.gz$|\\.fq\\.gz$|\\.fastq$|\\.fq$", "_val_1.fq.gz", pre.r1.base )
		full.trimmed.r1 <- paste(PARAMS[[ "FASTQ_TRIMMED_DIR" ]], trimmed.r1, sep="" )

		trimmed.r2 <- gsub("\\.fastq\\.gz$|\\.fq\\.gz$|\\.fastq$|\\.fq$", "_val_2.fq.gz", pre.r2.base )
		full.trimmed.r2 <- paste(PARAMS[[ "FASTQ_TRIMMED_DIR" ]], trimmed.r2, sep="" )

		}else
		{

		##it only validates if using paried end data
		trimmed.r1 <- gsub("\\.fastq\\.gz$|\\.fq\\.gz$|\\.fastq$|\\.fq$", "_trimmed.fq.gz", pre.r1.base )
		full.trimmed.r1 <- paste(PARAMS[[ "FASTQ_TRIMMED_DIR" ]], trimmed.r1, sep="" )

		full.trimmed.r2 <- "SINGLE_END"
		}

		
		if( !SUBMIT.FROM.ENGINE && !file.exists(full.trimmed.r1) )##this won't ever exist at submission time.
		{
		msg <- paste("Out of Coffee error in SV_ENGINE.R after trimming \"", this.lims.id, "\":\n\t",
					"Failed to correctly infer the file name for R1-trimmed.\n",
					"Tried:\n\t\"", full.trimmed.r1, "\"\n", sep="")
		stop(msg)
		}

		if( !SUBMIT.FROM.ENGINE && PARAMS[["FASTQ_2"]] != "SINGLE_END" && !file.exists(full.trimmed.r2) )
		{
		msg <- paste("Out of Coffee error in SV_ENGINE.R after trimming \"", this.lims.id, "\":\n\t",
					"Failed to correctly infer the file name for R2-trimmed.\n",
					"Tried:\n\t\"", full.trimmed.r2, "\"\n", sep="")
		stop(msg)
		}

	PARAMS[["TRIMMED_FASTQ_1"]] <- full.trimmed.r1
	PARAMS[["TRIMMED_FASTQ_2"]] <- full.trimmed.r2



	###########################
	###			###
	###	FASTQC_TRIMMED	###############################################################################################
	###			###
	###########################


fastqc.1.trimmed.done <- paste(PARAMS[["FASTQC_TRIMMED_DIR"]], ".", this.lims.id, "_R1.fastqc_trimmed.done", sep="")
fastqc.2.trimmed.done <- paste(PARAMS[["FASTQC_TRIMMED_DIR"]], ".", this.lims.id, "_R2.fastqc_trimmed.done", sep="")

	if( !file.exists(fastqc.1.trimmed.done) || FORCE.PROCESS.DOWNSTREAM )
	{
	##no need to set FORCE.PROCESS.DOWNSTREAM here as these results don't impact anything downstream

	this.input.fastq.file <- PARAMS[["TRIMMED_FASTQ_1"]]
	this.output.directory <- PARAMS[["FASTQC_TRIMMED_DIR"]]
	
	these.script.params <- paste(
					"\"FASTQC=", BINARY.CALLS[["FASTQC"]], "\" ",##defined in the config file
					"\"INPUT.FASTQ.FILE=", this.input.fastq.file, "\" ", 
					"\"OUTPUT.DIRECTORY=", this.output.directory, "\" ",
					"\"COMPLETION.FILE.NAME=", fastqc.1.trimmed.done, "\" ",
					sep="")


	this.job.name <- paste( "tFQCr1_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_r1_fqc_log.txt", sep="")

	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "FASTQC.R" ]], " ", these.script.params, sep="")


		if( SUBMIT.FROM.ENGINE )
		{

		deps <- NULL
			if( length(DEPENDENCIES) > 0 )
			{
			deps <- DEPENDENCIES
			}

		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)




		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (trimmed) for Read 1 with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()


		S <- system(sub.cmd, intern=TRUE)

		this.id <- get.scheduler.id(S) 
 			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}

		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)

		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)		##adding this to the dependencies list to avoid hitting 
									##the same fastq file with QC _and_ alignment at the same time

		msg <- paste("Job \"", this.id, "\" submitted:\n\t\"", S, "\"\n", sep="")
		cat(msg)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (trimmed) for Read 1 directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}


	}else
	{
	msg <- paste("\tFASTQC (trimmed) for Read 1 already done: skipping!\n")
	cat(msg)
	}

	if( !file.exists(fastqc.2.trimmed.done) || FORCE.PROCESS.DOWNSTREAM )
	{
	##no need to set FORCE.PROCESS.DOWNSTREAM here as these results don't impact anything downstream

	this.input.fastq.file <- PARAMS[["TRIMMED_FASTQ_2"]]
	this.output.directory <- PARAMS[["FASTQC_TRIMMED_DIR"]]
	
	these.script.params <- paste(
					"\"FASTQC=", BINARY.CALLS[["FASTQC"]], "\" ",##defined in the config file
					"\"INPUT.FASTQ.FILE=", this.input.fastq.file, "\" ", 
					"\"OUTPUT.DIRECTORY=", this.output.directory, "\" ",
					"\"COMPLETION.FILE.NAME=", fastqc.2.trimmed.done, "\" ",
					sep="")


	this.job.name <- paste( "tFQCr2_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_r2_fqc_log.txt", sep="")

	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "FASTQC.R" ]], " ", these.script.params, sep="")


		if( SUBMIT.FROM.ENGINE )
		{

		deps <- NULL
			if( length(DEPENDENCIES) > 0 )
			{
			deps <- DEPENDENCIES
			}

		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)



		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (trimmed) for Read 2 with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()
		S <- system(sub.cmd, intern=TRUE)

		this.id <- get.scheduler.id(S) 
 			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}
 
		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)


		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)		##adding this to the dependencies list to avoid hitting 
									##the same fastq file with QC _and_ alignment at the same time

		msg <- paste("Job \"", this.id, "\" submitted:\n\t\"", S, "\"\n", sep="")
		cat(msg)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling FastQC (trimmed) for Read 2 directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}


	}else
	{
	msg <- paste("\tFASTQC (trimmed) for Read 2 already done: skipping!\n")
	cat(msg)
	}

	###################
	###		###
	###	ALIGN	#########################################################################################################
	###		###
	###################
cat(">>>>>", date(), "\tSetting up alignment call on trimmed fastq files.\n", sep="")
flush.console()

aln.done <- paste(PARAMS[["BAM_SORTED_DIR"]], ".", this.lims.id, "_alignments.done", sep="")

this.sam <- paste(PARAMS[["SAM_DIR"]], this.lims.id, ".sam", sep="")
this.bam <- paste(PARAMS[["BAM_DIR"]], this.lims.id, ".bam", sep="")
this.sorted <- paste(PARAMS[["BAM_SORTED_DIR"]], this.lims.id, "_sorted.bam", sep="")
this.metric.file <- paste(PARAMS[["BAM_SORTED_DIR"]], this.lims.id, "_sorted_metrics.txt", sep="")
this.detailed.metric.file <- paste(PARAMS[["BAM_SORTED_DIR"]], this.lims.id, "_sorted_detailed_metrics.txt", sep="")
PARAMS[["SORTED_BAM"]] <- this.sorted
PARAMS[["SORTED_METRICS"]] <- this.metric.file

	if( !file.exists(aln.done) || FORCE.PROCESS.DOWNSTREAM)
	{
	FORCE.PROCESS.DOWNSTREAM <- TRUE

	this.output.directory <- PARAMS[["BAM_SORTED_DIR"]]	##one could argue that the log should be written somewhere else

	these.params <- paste(
		"\"ALN=", BINARY.CALLS[["ALN"]], "\" ", 
		"\"SAMTOOLS=", BINARY.CALLS[["SAMTOOLS"]], "\" ",
		"NTHREADS=", PARAMS[["MAX_THREADS"]], " ",
		"REFERENCE=", PARAMS[["REFERENCE_PATH"]], " ",
		"FASTQ_READ_1=", PARAMS[["TRIMMED_FASTQ_1"]], " ", 
		"FASTQ_READ_2=", PARAMS[["TRIMMED_FASTQ_2"]], " ",##set="SINGLE_END" if... single end data
		"OUTPUT_SAM_FILE=", this.sam, " ",
		"OUTPUT_BAM_FILE=", this.bam, " ", 
		"OUTPUT_SORTED_FILE=", this.sorted, " ", 
		"OUTPUT_METRICS_FILE=", this.metric.file, " ",
		"OUTPUT_DETAILED_METRICS_FILE=", this.detailed.metric.file, " ", 
		"COMPLETION.FILE.NAME=", aln.done,##file to be created to indicate _successful_ completion of all the jobs.
		sep="")


	this.job.name <- paste( "Aln_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_alignment_log.txt", sep="")
	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[["ALIGN.R"]], " ", 
 					these.params,
					sep="")


		if( SUBMIT.FROM.ENGINE )
		{

		deps <- NULL
			if( length(DEPENDENCIES) > 0 )
			{
			deps <- DEPENDENCIES

			}

		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)


		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling alignments with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()
		S <- system(sub.cmd, intern=TRUE)

		this.id <- get.scheduler.id(S) 
 			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}
  
		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)

		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling alignments directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}

	}else
	{
	msg <- paste("\tAlignments have already been done: skipping!\n")
	cat(msg)
	}


	###########################
	###			###
	###	QDNAseq	(no XY)	########################################################################################
	###			###
	###########################
cat(">>>>>", date(), "\tSetting up QDNAseq call on sorted BAM files.\n", sep="")
flush.console()

qdna.done <- paste(PARAMS[["QDNASEQ_DIR"]], ".", this.lims.id, "_qdnaseq.done", sep="")



	if( 	!file.exists( qdna.done ) || FORCE.PROCESS.DOWNSTREAM )
	{
	FORCE.PROCESS.DOWNSTREAM <- TRUE
	this.annotated.bins.file <- PARAMS[["ANNOTATED_BINS_FILE"]]
	these.filter.chromosomes <- "NA"				##how would we pass these at the commandline/in the design file?


	this.output.directory <- PARAMS[["QDNASEQ_DIR"]]

	out.pdf <- paste(PARAMS[["QDNASEQ_DIR"]], this.lims.id, "_qdnaseq-copy-number.pdf", sep="")
	out.rdat <- paste(PARAMS[["QDNASEQ_DIR"]], this.lims.id, "_qdnaseq-copy-number.RDat", sep="")

	these.params <- paste(
		"CONFIG.R_PATH=",  PARAMS[["CONFIG.R_PATH"]], " ",
		"SORTED_BAM=", PARAMS[["SORTED_BAM"]], " ",
		"SAMPLE_NAME=",this.lims.id, " ",
		"SAMPLE_DESCRIPTION=", this.sample.name.raw, " ",
		"BASE_BINS_FILE=", this.annotated.bins.file, " ",
		"ANNO_BED_FILE=", PARAMS[["QDNASEQ_BED_FILE"]], " ",
		"FILTER_CHROMOSOMES=", these.filter.chromosomes, " ",
		"OUTPUT_PDF=", out.pdf, " ",
		"OUTPUT_RDAT=", out.rdat, " ",

		"INCLUDE_SEX_CHROMOSOMES=FALSE", " ",
		"COMPLETION.FILE.NAME=", qdna.done, " ",
		sep="")



	this.job.name <- paste( "QDNA_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_qdna_log.txt", sep="")

	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "QDNASEQ.R" ]], " ", 
 					these.params,
					sep="")

		if( SUBMIT.FROM.ENGINE )
		{

		deps <- NULL
			if( length(DEPENDENCIES) > 0 )
			{
			deps <- DEPENDENCIES
			}

		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)


		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling QDNA (no XY) with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()
		S <- system(sub.cmd, intern=TRUE)

		this.id <- get.scheduler.id(S) 
 			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}

		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)
 
		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)##do we want to wait for this to finish before starting the run to include XY? 
								##could make it optional at some point but we'll inforce the wait for now.

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling QDNASeq (no XY) directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}

	}else
	{
	msg <- paste("\tQDNAseq has already been run: skipping!\n")
	cat(msg)
	}

	###################################
	###				###
	###	XY_QDNAseq		############################################################################################
	###				###
	###################################
cat(">>>>>", date(), "\tSetting up QDNAseq (XY) call on sorted BAM files.\n", sep="")
flush.console()

xy.qdna.done <- paste(PARAMS[["XY_QDNASEQ_DIR"]], ".", this.lims.id, "_xy.qdnaseq.done", sep="")

	if( !file.exists( xy.qdna.done ) || 	FORCE.PROCESS.DOWNSTREAM )
	{
	FORCE.PROCESS.DOWNSTREAM <- TRUE

	this.annotated.bins.file <- PARAMS[["ANNOTATED_BINS_FILE"]]
	these.filter.chromosomes <- "NA"

	this.output.directory <- PARAMS[["XY_QDNASEQ_DIR"]]

	out.pdf <- paste(PARAMS[["XY_QDNASEQ_DIR"]], this.lims.id, "_xy-qdnaseq-copy-number.pdf", sep="")
	out.rdat <- paste(PARAMS[["XY_QDNASEQ_DIR"]], this.lims.id, "_xy-qdnaseq-copy-number.RDat", sep="")

	these.params <- paste(
		"CONFIG.R_PATH=",  PARAMS[["CONFIG.R_PATH"]], " ",
		"SORTED_BAM=", PARAMS[["SORTED_BAM"]], " ",
		"SAMPLE_NAME=",this.lims.id, " ",
		"SAMPLE_DESCRIPTION=", this.sample.name.raw, " ",
		"BASE_BINS_FILE=", this.annotated.bins.file, " ",
		"ANNO_BED_FILE=", PARAMS[["QDNASEQ_BED_FILE"]], " ",
		"FILTER_CHROMOSOMES=", these.filter.chromosomes, " ",
		"OUTPUT_PDF=", out.pdf, " ",
		"OUTPUT_RDAT=", out.rdat, " ",

		"INCLUDE_SEX_CHROMOSOMES=TRUE", " ",
		"COMPLETION.FILE.NAME=", xy.qdna.done, " ",
		sep="")

	this.job.name <- paste( "xyQDNA_", this.lims.id, sep="")
	this.log.file <- paste(this.output.directory, this.lims.id, "_xyqdna_log.txt", sep="")

	this.job.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[[ "QDNASEQ.R" ]], " ", 
 					these.params,
					sep="")

		if( SUBMIT.FROM.ENGINE )
		{

		deps <- NULL
			if( length(DEPENDENCIES) > 0 )
			{
			deps <- DEPENDENCIES
			}

		sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)

		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling QDNA (with XY) with\n\t\"", sub.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()
		S <- system(sub.cmd, intern=TRUE)

		this.id <- get.scheduler.id(S) 
 			if( is.na(this.id) )
			{
			msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
			stop(msg)
			}
 
		msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
		cat(msg)

		add <- c( this.lims.id, this.job.name, this.id, "")
		SUBMISSION.IDS <- rbind(SUBMISSION.IDS, add)
		DEPENDENCIES <- c( DEPENDENCIES, this.id)

		}else
		{
		msg <- paste(">>>>>", date(), "\tSV_ENGINE calling QDNASeq (with XY) directly with\n\"", this.job.cmd, "\"\n", sep="")
		cat(msg)
		flush.console()

		S <- system(this.job.cmd, intern=FALSE)
			if( S != 0 )
			{
			stop("Job failed.\n")
			}
		}

	}else
	{
	msg <- paste("\tXY_QDNAseq has already been run: skipping!\n")
	cat(msg)
	}


	###################################################
	###						###
	###	Individual Sample Processing Completed	############################################
	###						###
	###################################################



##
##If we've been submitting stuff then we need to wait for the last job to return successfully before we return.
##
	if( SUBMIT.FROM.ENGINE )
	{

	##we'll submit a job to touch a completion file which depends on the current DEPENDENCIES 
	##	 and just sleep until it exists.

	
	##Once submitted, we could check the queue for the presence of its job id, sleep through intervals and once we find the 
	##	id to be missing check for the existence of the completion file


	this.job.cmd <- paste("touch ", PARAMS[["COMPLETION_FILE_NAME"]], sep="")
	this.job.name <- paste( "End_", this.lims.id, sep="")
	this.log.file <- paste(PARAMS[["RESULTS_DIR"]], this.lims.id, "_Completed.txt", sep="")

		
	deps <- NULL
		if( length(DEPENDENCIES) > 0 )
		{
		deps <- DEPENDENCIES
		}

	sub.cmd <- build.scheduler.submission(
					scheduler.options = PARAMS[["SCHEDULER_OPTIONS"]],
					dependencies = deps,
					job.name =this.job.name,
					log.file = this.log.file,
					job.cmd = this.job.cmd
					)

	msg <- paste(">>>>>", date(), "\tSV_ENGINE calling completion watcher with\n\t\"", sub.cmd, "\"\n", sep="")
	cat(msg)
	flush.console()

	S <- system(sub.cmd, intern=TRUE)

	this.id <- get.scheduler.id(S) 
		if( is.na(this.id) )
		{
		msg <- paste("Error submitting job for \"", this.lims.id, "\"\n\t\"", S, "\"\n", sep="")
		stop(msg)
		}
 	msg <-paste("Job submitted with id \"", this.id, "\"\n", sep="")
	cat(msg)

		while(!is.na( check.queue.for.id( this.id ) ) )
		{
		Sys.sleep(10*60)
		}
	

		if( file.exists(PARAMS[["COMPLETION_FILE_NAME"]]) )
		{

		cat("All finished for this sample @ ", date(), "\n", sep="")

		}else
		{
		msg <- paste("Detected a step fail.\n")
		stop(msg)
		}

	}else
	{

	
	cmd <- paste("touch ", PARAMS[["COMPLETION_FILE_NAME"]], sep="")
	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Ironically, creation of overall successful completion marker failed.\n\t\"", 
			PARAMS[["COMPLETION_FILE_NAME"]], "\"\n", sep="")
		stop(msg)
		}

	cat("All finished for this sample @ ", date(), "\n", sep="")

	}

setwd(init.dir)##this is obviously only useful if debugging...