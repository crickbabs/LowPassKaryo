
		###
		##THIS FILE WILL SELF-TEST IF SOURCED DIRECTLY FROM WITHIN R##
		##############################################################




			###################################################################################################################
			###														###
			###		PLEASE CUSTOMISE THE PARAMETERS BELOW FOR YOUR LOCAL SYSTEM PRIOR TO RUNNING THE PIPELINE	###
			###														###
			###################################################################################################################

			###################################################################
			###								###
			###	Primary Site-Specific Config & Default Parameters File	###	
			###								###
			###################################################################


CONTACT.ADDRESS <- ##this email address will be included at the top of the report file as a first point of contact for queries.

##central source directory
SOURCE.DIR <- ##absolute path to the directory containing the source files - by default it will be the directory containing this file


##need to document the script to generate this (or at least the file format)
##Absolute path to the genome lookup file -- see SetUp for details of the appropriate format for this file
GENOMES.LOOKUP.FILE <- ##



##
##The system specific call for each piece of software: basically, include here everything you would need to enter at the command line to invoke the tool, including loading 'modules' as appropriate.
##
BINARY.CALLS <- list(
			"RSCRIPT" = ##what would one type at the command line to start Rscript
##e.g. under a module system one might need to set
##	"RCRIPT" = "module purge ; module use /tools/easybuild/modules/all ;  module load R/3.5.1-foss-2016b ; Rscript ",
			"FASTQC" = ##Quality Control tool
			"QCTRIM" = ##the quality/adapter trimming tool. Existing code will work with trimmomatic.
			"ALN" = ##the alignment tool of choice - the existing code is set up to use bwa
			"SAMTOOLS" = ##command to run samtools
			)


BINARY.NAMES <- list(##what would we like to call the binary in the report/descriptions?
		"RSCRIPT" = "Rscript",
		"FASTQC" = "FastQC",
		"QCTRIM" = "Trimgalore",
		"ALN" = "bwa",
		"SAMTOOLS" = "samtools"
		)


VERSION.COMMANDS <- list(##how do we get the binary to report its version?
			"RSCRIPT" = " --version 2>&1 ",
			"FASTQC" = " --version ",
			"QCTRIM" = " --version ",
			"ALN" = " 2>&1 | grep Version" ,
			"SAMTOOLS" = " 2>&1 | grep Version"
			)

WRAPPED.FILES <- list(##wrapper scripts to run each step of the pipeline for a specific piece of software
			"SV_ENGINE.R" = "SV_ENGINE.R",
			"FASTQC.R" = "WRAPPED_do_fastqc.R",
			"QCTRIM.R" = "WRAPPED_do_trimgalore.R",
			"ALIGN.R" = "WRAPPED_do_bwa.R",
			"QDNASEQ.R" = "WRAPPED_do_qdnaseq.R",
			"PREP_REPORT.R" = "WRAPPED_do_build_project_summary.R",
			"CREATE_REPORT.R" = "WRAPPED_do_build_html_report.R"
			)


LOCAL.R2.REGEXP <- "_R2_001\\.fastq\\.gz$"	##what is the naming convention for your read 2 files?
						##this will be used to infer whether you have a single or paired end run
						
LOCAL.R1.REGEXP <- gsub("_R2_", "_R1_", LOCAL.R2.REGEXP)



##individual file names -- NOT absolute paths but relative to SOURCE.DIR
FILES.TO.BE.SOURCED <- list(
			UTILITY.FUNCTIONS.FILE = "utility_functions.R",
			BUILD.DESIGN.TEMPLATE.FILE = "build_design_template.R",
			BUILD.SCHEDULER.SUBMISSIONS.FILE = "scheduler_submission_functions.R"
			)



		###########################################################################################################################
		###															###
		###	The following default parameters do not _need_ to be altered for your location but you may wish to do so	###
		###															###
		###########################################################################################################################



DEFAULT.ANNO.BED.FILE <- NA
DEFAULT.ANNO.COLOUR <- "lightblue"
DEFAULT.MAX.THREADS <- 8
COMPLETION.MARKER.FILE.NAME <- ".all.done"	##this will appear in the results directory after a successful run

MINIMAL.COVERAGE.FLAG.THRESHOLD <- 50


##The following are used only in the generation of the template design file
DEFAULT.PROJECT.NAME <- "LowPassKaryo_Project_1"
DEFAULT.FASTQ.SOURCE.DIR <- "/path/to/my/example/fastq/files/"
DEFAULT.OUTPUT.ROOT.DIR <- "/path/to/my/first/test/run"



DEBUG <- FALSE
VERBOSE <- FALSE
SUBMIT.FROM.ENGINE <- FALSE	##option to submit jobs from the engine job - allows parallel processing of paired read files 
				##for FastQC and possible expansion downstream but if the cluster is busy then this just 
				##makes it more likely that jobs will get stuck in a queue. 
				##Running with this option set to TRUE is only lightly tested and is NOT recommended. 
				## Default (FALSE) is off.
				#######################################################################################





			###########################################################
			###########################################################
			###########################################################
			###							###
			###							###
###########################		DO NOT EDIT BELOW THIS LINE		################################
			###							###
			###							###
			###########################################################
			###########################################################
			###########################################################

VALID.COMMAND.LINE.OPTIONS <- list(
				"--help" = list(
						"DEFAULT" = NA,
						"Description" = "Print this help message."
						),
                		"--available-genomes" = list(
                                				"DEFAULT" = NA,
                                				"Description" = "Print a list of currently available genomes along with their \"Tag\" IDs."
                                				),
                		"--template-design" = list(
 	                               			"DEFAULT" = "design_template_file.txt",
        			                        "Description" = paste("Generate a file in the format expected for a design file:\n",
                                                        		"\t\tThe file will be written to the name passed with the argument.",sep="")
                                			),
                		"--design" = list(
                        		        "DEFAULT" = "design_file_name.txt",
                                		"Description" = "Name of the design file to be used for this run."
                                		)#,
				##this is going to be in the design file.
				#"--scheduler-options" = list(
				#		"DEFAULT" = "",
				#		"Description" = "Options to be passed to the scheduler:\n\t\tIf used, this will replace _all_ of the default option set."
				#		)
                		)



	############
	##PATH.CAP##
	############
########################################################################################################################
path.cap <- function(
        input.path
        )
{

N <- nchar(input.path)

slash <- .Platform$file.sep

        if( substr(x = input.path, start = N, stop=N) != slash )
        {
        input.path <- paste(input.path, slash, sep="")
        }

return(input.path)
}##end PATH.CAP
########################################################################################################################



SOURCE.DIR <- path.cap(SOURCE.DIR)
ERRORS <- vector()


		####################
		##CHECK.FILE.PATHS##
		####################
########################################################################################################################
check.file.paths <- function(X)
{
	##
	##Convert to absolute paths
	##
	for(i in seq_along(X))
	{
	this.name <- names(X)[i]
	this.file <- X[[ this.name ]]
	full.file <- paste(SOURCE.DIR, this.file, sep="")
		if( !file.exists(full.file) )
		{
		ERRORS <- c( ERRORS, paste("Provided \"", this.name, 
						"\" file could not be located in the SOURCE.DIR: looked for\n\t\"", 
						this.file, "\" in \n\t\"", SOURCE.DIR, "\"\n", sep=""))
		}

	X[[ this.name ]] <- full.file
	}



	if( length(ERRORS) > 0 )
	{
	stop(ERRORS)
	}

return(X)
}##end CHECK.FILE.PATHS
########################################################################################################################



	if( VERBOSE || DEBUG )
	{
	cat("Performing source file sanity check.\n")
	}
FILES.TO.BE.SOURCED <- check.file.paths(X=FILES.TO.BE.SOURCED)

	if( VERBOSE || DEBUG )
	{
	cat("\tDone\n")
	}

	for(i in seq_along(FILES.TO.BE.SOURCED) )
	{
	source(FILES.TO.BE.SOURCED[[i]])
	}

	if( !exists("build.scheduler.submission") )
	{
	msg <- paste("Please ensure that the function \"build.scheduler.submission\" is defined in the \"BUILD.SCHEDULER.SUBMISSIONS.FILE\" file:\n\t\"",
				FILES.TO.BE.SOURCED[["BUILD.SCHEDULER.SUBMISSIONS.FILE"]], "\"\n", sep="")
	stop(msg)
	}
	if( !exists( "DEFAULT.SCHEDULER.OPTIONS" ) )
	{
	msg <- paste("Please ensure that the variable \"DEFAULT.SCHEDULER.OPTIONS\" is defined in the \"BUILD.SCHEDULER.SUBMISSIONS.FILE\" file:\n\t\"",
				FILES.TO.BE.SOURCED[["BUILD.SCHEDULER.SUBMISSIONS.FILE"]], "\"\n", sep="")
	stop(msg)
	}

	if( VERBOSE || DEBUG )
	{
	cat("Performing wrapped file sanity check.\n")
	}

WRAPPED.FILES <- check.file.paths(X=WRAPPED.FILES)
	if( VERBOSE || DEBUG )
	{
	cat("\tDone\n")
	}



w <- which(!names(BINARY.CALLS) %in% names(BINARY.NAMES))
	if( length(w) > 0 )
	{
	msg <- paste("The following BINARY.CALLS entries had no corresponding BINARY.NAMES entry:\n\t", 
			paste(names(BINARY.CALLS)[w], collapse="\n\t"), "\n", sep="")
	stop(msg)
	}

w <- which(!names(BINARY.CALLS) %in% names(VERSION.COMMANDS))
	if( length(w) > 0 )
	{
	msg <- paste("The following BINARY.CALLS entries had no corresponding VERSION.COMMANDS entry:\n\t", 
			paste(names(BINARY.CALLS)[w], collapse="\n\t"), "\n", sep="")
	stop(msg)
	}




