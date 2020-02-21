


		###########################
		###			###
		###	DO_BWA		###
		###			###
		###########################

	###########################################################################
	##This script is intended to be automatically invoked by a wrapper script##
	###########################################################################



my.args <- commandArgs( trailingOnly = TRUE )

EXPECTED.PARAMS <- c(
		"ALN",		##system specific method for invoking FASTQC (e.g. including module load statements on camp)
		"SAMTOOLS",
		"NTHREADS",
		"REFERENCE",
		"FASTQ_READ_1",
		"FASTQ_READ_2",	##set="SINGLE_END" if... single end data
		"OUTPUT_SAM_FILE",
		"OUTPUT_BAM_FILE",
		"OUTPUT_SORTED_FILE",
		"OUTPUT_METRICS_FILE",
		"OUTPUT_DETAILED_METRICS_FILE",
		"COMPLETION.FILE.NAME"	##file to be created to indicate _successful_ completion of all the jobs.
		)

VALID.OTHER.PARAMETERS <- vector()




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

MSG <- paste("\n\nWRAPPED_do_bwa.R invoked on \"", this.host, "\"\n\t@", d, "\n\n", sep="")
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


	if( !file.exists( PARAMS[["FASTQ_READ_1"]] ) )
	{
	msg <- paste("Unable to locate provided R1 fastq file:\n\t\"", PARAMS[["FASTQ_READ_1"]], "\"\n", sep="")
	stop(msg)
	}

	if( PARAMS[[ "FASTQ_READ_2" ]] != "SINGLE_END" && !file.exists( PARAMS[["FASTQ_READ_2"]] ) )
	{
	msg <- paste("Unable to locate provided R2 fastq file:\n\t\"", PARAMS[["FASTQ_READ_2"]], "\"\n", sep="")
	stop(msg)
	}

PARAMS[["BWA"]] <- PARAMS[["ALN"]]

	###########################
	###			###
	###	FASTQ2SAM	###
	###			###
	###########################


this.sam <- PARAMS[["OUTPUT_SAM_FILE"]]
this.sam <- unlist(strsplit(this.sam, split=.Platform$file.sep))
L <- length(this.sam)
sam.done <- paste(c(this.sam[1:(L-1)], paste(".", this.sam[L], ".done", sep="")), collapse=.Platform$file.sep)



	if( !file.exists(sam.done) )
	{

	cat("Performing alignment:\n")
	flush.console()
	RG <- paste("\'@RG\\tID:",  PARAMS[["FASTQ_READ_1"]], 
						"\\tDS:Short",##this is required by CAVEMAN (and probably PICARD by now)
						"\\tPL:ILLUMINA", ##Hard coding this because GATK doesn't recognise IonTorrent so it's pointless to change it.
						"\\tSM:", PARAMS[["FASTQ_READ_1"]], ##for a single run pool & sample can be identical -- merging will re-header anyway
						"\'", sep="")


	bwa.cmd <- paste(
			PARAMS[["BWA"]], " mem ", 
				" -R ", RG, 
				" -t ", PARAMS[[ "NTHREADS" ]],
				" -M ", ##mark secondary hits so that Picard doesn't throw a hissy fit.
			##
			##Keep the rest at default for now
			##
			#####Algorithm options:

			       " -k 19 ", 	#INT        minimum seed length [19]
			       " -w 100 ", 	#INT        band width for banded alignment [100]
			       " -d 100 ", 	#INT        off-diagonal X-dropoff [100]
			       " -r 1.5 ", 	#FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
			       " -y 20 ", 	#INT        seed occurrence for the 3rd round seeding [20]
			       " -c 500 ", 	#INT        skip seeds with more than INT occurrences [500]
			       " -D 0.5 ", 	#FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
			       " -W 0 ", 	#INT        discard a chain if seeded bases shorter than INT [0]
			       " -m 50 ", 	#INT        perform at most INT rounds of mate rescues for each read [50]
			       ##" -S ",          #skip mate rescue
			       ##" -P  ",         # skip pairing; mate rescue performed unless -S also in use
			#####Scoring options:
			       " -A 1 ",	#INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
			       " -B 4 ", 	#INT        penalty for a mismatch [4]
			       " -O 6,6 ", 	#INT[,INT]  gap open penalties for deletions and insertions [6,6]
			       " -E 1,1 ", 	#INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
			       " -L 5,5 ",	#INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
			       " -U 17 ",	#INT        penalty for an unpaired read pair [17]
	
		     #  -x STR        read type. Setting -x changes multiple parameters unless overriden [null]
		     #                pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
		     #                ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
		     #                intractg: -B9 -O16 -L5  (intra-species contigs to ref)

				" ", PARAMS[[ "REFERENCE" ]],
				" ", PARAMS[[ "FASTQ_READ_1" ]], 
				sep="")

		if( PARAMS[[ "FASTQ_READ_2" ]] != "SINGLE_END" )
		{
		bwa.cmd <- paste(bwa.cmd, " ", PARAMS[[ "FASTQ_READ_2" ]], " ", sep="")
		}


	bwa.cmd <- paste(bwa.cmd, " > ", PARAMS[["OUTPUT_SAM_FILE" ]], sep="")

	cat("Calling bwa with\n\t\"", bwa.cmd, "\"\n", sep="")
	S <- system(bwa.cmd, intern=FALSE)
		if( S != 0 )
		{
		stop("bwa alignment failed.\n")
		}
	cmd <- paste("touch ", sam.done, sep="")
	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Ironically, creation of successful completion marker failed.\n\t\"", sam.done, "\"\n", sep="")
		stop(msg)	
		}

	}


	###################
	###		###
	###	SAM2BAM	###
	###		###
	###################


this.bam <- PARAMS[["OUTPUT_BAM_FILE"]]
this.bam <- unlist(strsplit(this.bam, split=.Platform$file.sep))
L <- length(this.bam)
bam.done <- paste(c(this.bam[1:(L-1)], paste(".", this.bam[L], ".done", sep="")), collapse=.Platform$file.sep)


	if( !file.exists(bam.done) )
	{
	input.sam <- PARAMS[["OUTPUT_SAM_FILE"]]
	output.bam <- PARAMS[["OUTPUT_BAM_FILE"]]

	bam.cmd <- paste(PARAMS[["SAMTOOLS"]], " view ", 
					" -S -b ",
					input.sam, " > ", output.bam,
					sep="")

	cat("Converting from Sam to Bam:\n\t\"", bam.cmd, "\"\n", sep="")
	S <- system(bam.cmd,intern=FALSE)
		if( S != 0 )
		{
		stop("Conversion failed.\n")
		}

	cmd <- paste("touch ", bam.done, sep="")
	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Ironically, creation of successful completion marker failed.\n\t\"", bam.done, "\"\n", sep="")
		stop(msg)	
		}
	}


	###################
	###		###
	###	SORTBAM	###
	###		###
	###################


this.sorted <- PARAMS[["OUTPUT_SORTED_FILE"]]
this.sorted <- unlist(strsplit(this.sorted, split=.Platform$file.sep))
L <- length(this.sorted)

sorted.done <- paste(c(this.sorted[1:(L-1)], paste(".", this.sorted[L], ".done", sep="")), collapse=.Platform$file.sep)

	if( !file.exists( sorted.done ) )
	{

	input.bam <- PARAMS[["OUTPUT_BAM_FILE"]]
	output.sorted <- PARAMS[["OUTPUT_SORTED_FILE"]]


	sort.cmd <- paste(PARAMS[["SAMTOOLS"]], " sort ",
					input.bam, 
					" ", gsub("\\.bam$", "", output.sorted), " ",
					##" > ", output.sorted,
					sep="")
					
	cat("Sorting bam file:\n\t\"", sort.cmd, "\"\n", sep="")
	flush.console()
	S <- system(sort.cmd, intern=FALSE)
		if( S != 0 )
		{
		stop("Sorting failed.\n")
		}
	
	cmd <- paste("touch ", sorted.done, sep="")
	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Ironically, creation of successful completion marker failed.\n\t\"", sorted.done, "\"\n", sep="")
		stop(msg)	
		}
	}


indexed.done <- paste(c(this.sorted[1:(L-1)], paste(".", this.sorted[L], ".bai.done", sep="")), collapse=.Platform$file.sep)
	if( !file.exists( indexed.done ) )
	{

	this.sorted <- PARAMS[["OUTPUT_SORTED_FILE"]]
	
	index.cmd <- paste(PARAMS[["SAMTOOLS"]], " index ", this.sorted, sep="")
	cat("Indexing sorted bam file:\n\t\"", index.cmd, "\"\n", sep="")
	flush.console()

	S <- system(index.cmd, intern=FALSE)
		if( S != 0 )
		{
		stop("Indexing failed.\n")

		}

	cmd <- paste("touch ", indexed.done, sep="")
	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Ironically, creation of successful completion marker failed.\n\t\"", indexed.done, "\"\n", sep="")
		stop(msg)	
		}

	}


	if( !file.exists( PARAMS[["OUTPUT_METRICS_FILE"]] ) )
	{

	this.sorted <- PARAMS[["OUTPUT_SORTED_FILE"]]

	idx.cmd <- paste(PARAMS[["SAMTOOLS"]], " idxstats ", this.sorted, sep="")
	cat("Extracting IDX Stats.\n\t\"", idx.cmd, "\"\n", sep="")
	flush.console()

	
	S <- system(idx.cmd, intern=TRUE)##we _do_ want to capture this output
	SS <- vecsplit(S, split="\t")

	HITS <- as.numeric(gsub(" +", "", SS[,3]))
	MISSES <- as.numeric(gsub(" +", "", SS[,4]))

	L <- length(HITS)

		if( SS[L,1] == "*" )
		{
		MISSES <- c(MISSES, HITS[L])
		HITS <- HITS[-L]
		}

	ALN <- sum(HITS)
	TOT <- ALN + sum(MISSES)

	CN <- c( "Sample", "Raw_Yield", "Tot_Aln", "Percent_Aln")
	REPORT <- matrix(
				c( this.sorted, TOT, ALN, round(100*ALN/TOT, 2)), 
				nrow = 1
			)

	colnames(REPORT) <- CN

	write.table(REPORT, file=PARAMS[["OUTPUT_METRICS_FILE"]], sep="\t", row.names=FALSE, col.names=TRUE)

	}


	if( !file.exists( PARAMS[["OUTPUT_DETAILED_METRICS_FILE"]] ) )
	{
	this.sorted <- PARAMS[[ "OUTPUT_SORTED_FILE" ]]

	idx.cmd <- paste(PARAMS[["SAMTOOLS"]], " idxstats ", this.sorted, sep="")
	cat("Extracting Detailed IDX Stats.\n\t\"", idx.cmd, "\"\n", sep="")
	flush.console()

	
	S <- system(idx.cmd, intern=TRUE)##we _do_ want to capture this output
	SS <- vecsplit(S, split="\t")

	REPORT.CN <- c("Chromosome", "IDXStat_Reported_Mapped_Reads")
	REPORT <- cbind(SS[,1], gsub(" +", "", SS[,3]) )
	colnames(REPORT) <- REPORT.CN
	
	write.table(REPORT, file=PARAMS[["OUTPUT_DETAILED_METRICS_FILE"]], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
	}


##all done.


cmd <- paste("touch ", PARAMS[["COMPLETION.FILE.NAME"]], sep="")
S <- system(cmd, intern=FALSE)
	if( S != 0 )
	{
	msg <- paste("Ironically, creation of overall successful completion marker failed.\n\t\"", PARAMS[["COMPLETION.FILE.NAME"]], "\"\n", sep="")
	stop(msg)	
	}
cat("bwa done.\n")