



		###########################
		###			###
		###	DO_QDNASEQ	###
		###			###
		###########################

	###########################################################################
	##This script is intended to be automatically invoked by a wrapper script##
	###########################################################################


my.args <- commandArgs(trailingOnly = TRUE )


	if( FALSE )
	{
	#cat("Setting parameters for TESTING\n")
	#my.args <- c()
	}


EXPECTED.PARAMS <- c(
		"CONFIG.R_PATH",
		"SORTED_BAM",
		"SAMPLE_NAME",
		"SAMPLE_DESCRIPTION",
		"BASE_BINS_FILE",
		"ANNO_BED_FILE",
		"FILTER_CHROMOSOMES",
		"OUTPUT_PDF",
		"OUTPUT_RDAT",
		"INCLUDE_SEX_CHROMOSOMES",
		"COMPLETION.FILE.NAME"
		)

VALID.OTHER.PARAMS <- vector()



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

MSG <- paste("\n\nWRAPPED_do_qdnaseq.R invoked on \"", this.host, "\"\n\t@", d, "\n\n", sep="")
cat(MSG)

cat("Attempting to load required R libraries:\n")

require(QDNAseq)	##this _is_ installed for module load R/3.5.1-foss-2016b-BABS
require(Biobase)


cat("\n\nCalled with the following parameters:\n\t",
			paste(my.args, collapse="\n\t"),
			"\n\n\n", sep="")



##As long as these are working, we can parse the command line params


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



PARAMS <- list()

	for(i in seq_along(EXPECTED.PARAMS))
	{
	this.param <- EXPECTED.PARAMS[i]
	PARAMS[[ this.param ]] <- tmp[M[i],2]
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

source(PARAMS[["CONFIG.R_PATH"]])

###Import the base_bins_table

cat("\tImporting the bin annotation file:\n")
load(PARAMS[["BASE_BINS_FILE"]])##loads an object called annotated.bins & default.chromosomes.to.filter

cat("\tImporting & binning read counts.\n")
binned.read.counts <- binReadCounts(annotated.bins, bamfile=PARAMS[["SORTED_BAM"]],
			bamnames=PARAMS[["SAMPLE_DESCRIPTION"]]
			##if we pass this a bamnames vector it'll prevent the file names being used to infer the sample name.
			)


	if( !is.na(PARAMS[["FILTER_CHROMOSOMES"]]) && 
		PARAMS[["FILTER_CHROMOSOMES"]] != "NA" &&
		 PARAMS[["FILTER_CHROMOSOMES"]] != " ")
	{
	filter.chromosomes <- PARAMS[["FILTER_CHROMOSOMES"]]

	}else
	{
	filter.chromosomes <- default.chromosomes.to.filter
	}


cat("\tApplying Filter:\n")
cat("\t\tExcluded chromosomes = \n\t", paste(filter.chromosomes, collapse="\n\t"), "\n", sep="")
pre.filtered.binned.read.counts <- applyFilters(
				object = binned.read.counts,
				residual = FALSE,
				blacklist = FALSE,
				mappability = TRUE,
				bases = TRUE,
				chromosomes = filter.chromosomes[! filter.chromosomes %in% c("X", "Y")]
				)

not.sex.filters <- fData(pre.filtered.binned.read.counts)[, "use"]


filtered.binned.read.counts <- applyFilters(
				object = binned.read.counts,
				residual = FALSE,
				blacklist = FALSE,
				mappability = TRUE,
				bases = TRUE,
				chromosomes = filter.chromosomes
				)

w <- which(pData(filtered.binned.read.counts)[,"used.reads"] == 0 )
	if( length(w) > 0 )
	{
	msg <- paste("No reads left for \n\t", 
					paste(row.names(pData)[w], collapse="\n\t"),
					"\n", sep="")

	stop(msg)
	}


cat("\tCalling estimateCorrection:\n")
corrected.filtered.binned.read.counts <- estimateCorrection(filtered.binned.read.counts)


	if( PARAMS[[ "INCLUDE_SEX_CHROMOSOMES" ]] == "TRUE" )
	{
	##
	##Add the sex chromosomes back in
	#################################
	cat("Un-masking sex chromosomes.\n")
	fData(corrected.filtered.binned.read.counts)[, "use"] <- not.sex.filters
	}

cat("\tCalling correctBins:\n")
copy.numbers <- correctBins(corrected.filtered.binned.read.counts)


cat("\tCalling normalizeBins:\n")
normalised.copy.numbers <- normalizeBins(copy.numbers)

cat("\tCalling smoothOutlierBins:\n")
smoothed.normalised.copy.numbers <- smoothOutlierBins(normalised.copy.numbers)

cat("\tCalling segmentBins:\n")
segmented.smoothed.normalised.copy.numbers <- segmentBins( smoothed.normalised.copy.numbers, transformFun="sqrt")	##"log2")

cat("\tCalling normalizeSegmentedBins:\n")
final.copy.numbers <- normalizeSegmentedBins( segmented.smoothed.normalised.copy.numbers)


	if( !is.na(PARAMS[["ANNO_BED_FILE"]]) && PARAMS[["ANNO_BED_FILE"]] != "NA" )
	{
	cat("Importing annotation bed file:\n\t\"", PARAMS[["ANNO_BED_FILE"]], "\n", sep="")
	this.bed <- as.matrix(read.delim(PARAMS[["ANNO_BED_FILE"]], sep="\n", header=FALSE, comment =""))
		if( length(grep("^#", this.bed[1])) == 1 )
		{
		this.bed <- this.bed[-1]
		}

	this.bed <- vecsplit(this.bed, split="\t")
	
	}else
	{
	this.bed <- matrix(NA, nrow= 0, ncol = 1)
	}


R <- length(PARAMS[["SORTED_BAM"]])
this.pdf.name <- PARAMS[["OUTPUT_PDF"]]
	if( nchar(this.pdf.name) == 0 )
	{
	stop("Empty pdf name.\n")
	}
cat("\tPlotting to \"", PARAMS[["OUTPUT_PDF"]], "\"\n", sep="")
pdf(PARAMS[["OUTPUT_PDF"]], width = 2*7, height = R*7)

plot(final.copy.numbers)

	if( length(this.bed) > 0 )
	{
	cat("\tAdding annotation to segmented plot:\n")
		##should this function be defined here (which prevents its use in the sumamry file)
		##or in utility_functions.R where it currently lives btu seems disconnected...
	RN <- 	add.bed.anno.to.qdna.plot(
				qdna.object = final.copy.numbers,
				bed = this.bed
				)
		
	}



dev.off()



cat("\tSaving R objects.\n")
save(final.copy.numbers, file=PARAMS[["OUTPUT_RDAT"]])

cmd <- paste("touch ", PARAMS[["COMPLETION.FILE.NAME"]], sep="")
S <- system(cmd, intern=FALSE)
	if( S != 0 )
	{
	msg <- paste("Ironically, creation of overall successful completion marker failed.\n\t\"", 
		PARAMS[["COMPLETION.FILE.NAME"]], "\"\n", sep="")
	stop(msg)
	}


cat("Done!\n")
