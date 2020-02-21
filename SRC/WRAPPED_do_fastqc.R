


		###########################
		###			###
		###	DO_FASTQC	###
		###			###
		###########################

	###########################################################################
	##This script is intended to be automatically invoked by a wrapper script##
	###########################################################################



my.args <- commandArgs( trailingOnly = TRUE )

EXPECTED.PARAMS <- c(
		"FASTQC",		##system specific method for invoking FASTQC (e.g. including module load statements on camp)
		"INPUT.FASTQ.FILE",	## fastq file you wish to QC
		"OUTPUT.DIRECTORY",	##where should FastQC write its results?
		"COMPLETION.FILE.NAME"	##file to be created to indicate _successful_ completion of the job.
		)

VALID.OTHER.PARAMS <- list(###taken from the -help report for current version of FASTQC.
##do it as a list to make it possibel to set default sif _just_ passed the switch name and no argument (or to make a full explicit call)
    "--casava" =NULL,		##Files come from raw casava output. Files in the same sample
                    #group (differing only by the group number) will be analysed
                    #as a set rather than individually. Sequences with the filter
                    #flag set in the header will be excluded from the analysis.
                    #Files must have the same names given to them by casava
                    #(including being gzipped and ending with .gz) otherwise they
                    #won't be grouped together correctly.

    "--nano" =NULL,		##Files come from nanopore sequences and are in fast5 format. In
                    #this mode you can pass in directories to process and the program
                    #will take in all fast5 files within those directories and produce
                    #a single output file from the sequences found in all files.

    "--nofilter" =NULL,      ##If running with --casava then don't remove read flagged by
                    ##casava as poor quality when performing the QC analysis.

    "--extract" =NULL,       ##If set then the zipped output file will be uncompressed in
                    #the same directory after it has been created.  By default
                    #this option will be set if fastqc is run in non-interactive
                    #mode.

    "-j" =NULL, "--java" =NULL,       ##Provides the full path to the java binary you want to use to
                    #launch fastqc. If not supplied then java is assumed to be in
                    #your path.

    "--noextract" =NULL,     ##Do not uncompress the output file after creating it.  You
                    #should set this option if you do not wish to uncompress
                    #the output when running in non-interactive mode.

    "--nogroup" =NULL,       ##Disable grouping of bases for reads >50bp. All reports will
                    #show data for every base in the read.  WARNING: Using this
                    #option will cause fastqc to crash and burn if you use it on
                    #really long reads, and your plots may end up a ridiculous size.
                    #You have been warned!

    "--min_length" =NULL,    ##Sets an artificial lower limit on the length of the sequence
                    #to be shown in the report.  As long as you set this to a value
                    #greater or equal to your longest read length then this will be
                    #the sequence length used to create your read groups.  This can
                    #be useful for making directly comaparable statistics from
                    #datasets with somewhat variable read lengths.

    "-f" =NULL,  "--format" =NULL,     ##Bypasses the normal sequence file format detection and
                    #forces the program to use the specified format.  Valid
                    #formats are bam,sam,bam_mapped,sam_mapped and fastq

    "-t" =NULL,  "--threads" =NULL,    ##Specifies the number of files which can be processed
                    #simultaneously.  Each thread will be allocated 250MB of
                    #memory so you shouldn't run more threads than your
                    #available memory will cope with, and not more than
                    #6 threads on a 32 bit machine

    "-c" =NULL,              #Specifies a non-default file which contains the list of
    "--contaminants" =NULL,  #contaminants to screen overrepresented sequences against.
                    #The file must contain sets of named contaminants in the
                    #form name[tab]sequence.  Lines prefixed with a hash will
                    #be ignored.

    "-a" =NULL,              #Specifies a non-default file which contains the list of
    "--adapters" =NULL,      #adapter sequences which will be explicity searched against
                    #the library. The file must contain sets of named adapters
                    #in the form name[tab]sequence.  Lines prefixed with a hash
                    #will be ignored.

    "-l" =NULL,              #Specifies a non-default file which contains a set of criteria
    "--limits" =NULL,        #which will be used to determine the warn/error limits for the
                    #various modules.  This file can also be used to selectively
                    #remove some modules from the output all together.  The format
                    #needs to mirror the default limits.txt file found in the
                    #Configuration folder.

   "-k" =NULL, "--kmers" =NULL,       #Specifies the length of Kmer to look for in the Kmer content
                    #module. Specified Kmer length must be between 2 and 10. Default
                    #length is 7 if not specified.

   "-q" =NULL, "--quiet" =NULL,       #Supress all progress messages on stdout and only report errors.

   "-d" =NULL, "--dir" =NULL         #Selects a directory to be used for temporary files written when
                    #generating report images. Defaults to system temp directory if
                    #not specified.
	)


	############
	##VECSPLIT##
	############
########################################################################################################################
vecsplit <- function(X, split="_",
                        cols = NULL)
{
g <- match(split, X)##stops after first match found rather than searching for all
        if( is.na(g) == 0 )
        {
        stop("Error split term not found\n")
        }


        if( is.null(cols) )
        #infer from the first row...
        {
        x <- length(unlist(strsplit(X[1], split=split)))
        cols <- 1:x
        }


x <- strsplit(X, split=split)

X <- matrix( NA, nrow = length(X), ncol = length(cols) )

        for( i in seq_along(cols) )
        {
        X[,i] <- unlist(lapply(x, "[", cols[i]))
        }

return(X)
}##end VECSPLIT
########################################################################################################################


d <- date()
this.host <- system("hostname", intern = TRUE)

MSG <- paste("\n\ndo_fastqc.R invoked on \"", this.host, "\"\n\t@", d, "\n\n", sep="")
cat(MSG)

	##
	##Check arguments, which should be passed in the form <NAME_1>=<VALUE> <NAME_2>=<VALUE> etc
	##
	if( length(my.args) == 0 )
	{
	msg <- paste("This script should be invoked with the following (named) arguments:\n\t",
					paste(EXPECTED.PARAMS, collapse="\n\t"), "\n", sep="")
	stop(msg)
	}	


tmp <- vecsplit(my.args, split="=")
M <- match( EXPECTED.PARAMS, tmp[,1] )
m <- which(is.na(M))
	if( length(m) > 0 )
	{
	msg <- paste("The following command line parameter names were not detected:\n\t",
				paste(EXPECTED.PARAMS[m], collapse="\n\t"), "\n", sep="")
	stop(msg)
	}

##allow for the passing of other, optional parameters, which will be introduced, _as_is_ into the system call
other.params <- tmp[-M,,drop=FALSE]
	if( length(other.params) > 0 )
	{
	w <- which(! other.params[,1] %in% names(VALID.OTHER.PARAMS) )
		if( length(w) > 0 )
		{
		msg <- paste("The following optional parameters were not recognised as valid:\n\t",
						paste(other.params[w,1], collapse="\n\t"), "\n", sep="")

		stop(msg)
		}

	msg <- paste("Detected ", dim(other.params)[1], " additional parameters.\n", sep="")
	cat(msg)
	}

PARAMS <- list()

	for(i in seq_along(EXPECTED.PARAMS))
	{
	this.param <- EXPECTED.PARAMS[i]
	PARAMS[[ this.param ]] <- tmp[M[i],2]
	}


	if( !file.exists( PARAMS[["INPUT.FASTQ.FILE"]] ) )
	{
	msg <- paste("Unable to locate provided fastq file:\n\t\"", PARAMS[["INPUT.FASTQ.FILE"]], "\"\n", sep="")
	stop(msg)
	}


	if( !file.exists( PARAMS[[ "OUTPUT.DIRECTORY" ]]  ) )
	{##all directories to be created _prior_ to calling to avoid sitting in the queue only to fail because the path's wrong.

	msg <- paste("Provided output directory does not exist:\n\t\"", PARAMS[[ "OUTPUT.DIRECTORY"]], "\"\n", sep="")
	stop(msg)
	}




##explicitly report the version (and, at the same time, check that the invocation actually works.
ver.cmd <- paste(PARAMS[["FASTQC"]] , " --version ", sep="")
ver <- system(ver.cmd, intern = TRUE)
g <- grep("-bash: fastqc: command not found", ver)
	if( length(g) > 0 )
	{
	msg <- paste("Attempt to invoke FastQC with the following command failed:\n\t\"", ver.cmd, "\"\n", sep="")
	stop(msg)
	}
msg <- paste("Performing FastQC using:\n\n\"", ver, "\"\n", sep="")
cat(msg)


	##############################
	###Build the system command###
	##############################

fastqc.cmd <- paste(
			PARAMS[["FASTQC"]], " ",
			" --outdir ", PARAMS[["OUTPUT.DIRECTORY"]], " ", 
			sep="")

	for(i in seq_along(other.params[,1]))##even if empty, should still have a dimension
	{
	fastqc.cmd <- paste( fastqc.cmd, " ", other.params[i,1], " ", other.params[i, 2], " ", sep="")
	msg <- paste("Added optional argument \"", other.params[i,1], "\" = \"",other.params[i,2], "\"\n", sep="")
	cat(msg)
	}

fastqc.cmd <- paste(fastqc.cmd, " ", PARAMS[[ "INPUT.FASTQ.FILE" ]], sep="")


msg <- paste("\n\n\tCalling FastQC:\t", date(), "\n", sep="")
msg <- paste(msg, "\n\n\t==============:\n\n", sep="")
msg <- paste(msg, "cmd:\n\t\"", fastqc.cmd, "\"\n", sep="")

cat(msg)

S <- system(fastqc.cmd, intern = FALSE)
	if( S != 0 )
	{
	msg <- "FastQC call failed.\n"
	stop(msg)

	}else
	{
	cmd <- paste("touch ", PARAMS[["COMPLETION.FILE.NAME"]], sep="")
	S <- system(cmd, intern = FALSE )
		if( S != 0 )
		{
		msg <- paste("Ironically, creation of successful completion marker failed.\n\t\"", PARAMS[["COMPLETION.FILE.NAME"]], "\"\n", sep="")
		stop(msg)
		}
	msg <- paste("All finished @ \"", date(), "\"\n", sep="")
	cat(msg)
	}
