


		###
		###This script will gather together all of the output metrics for a given project
		###	and generate summary files.
		###

my.args <- commandArgs(trailingOnly = TRUE)


EXPECTED.PARAMS <- c(
		"CONFIG.R_PATH",	
		"PARAMS_LOOKUP_FILE",	##we should be able to get everythign from here rather than re-constructing everything
		"COMPLETION_FILE_NAME"
		)

	if(FALSE )
	{
	#cat("Setting parameters for TESTING.\n")
	#my.args <- c()
	}
VALID.OTHER.PARAMETERS <- vector()



d <- date()
this.host <- system("hostname", intern = TRUE)

MSG <- paste("\n\nWRAPPED_do_build_project_sumamry.R invoked on \"", this.host, "\"\n\t@", d, "\n\n", sep="")
cat(MSG)


cat("\n\nCalled with the following parameters:\n\t", paste(my.args, collapse="\n\t"), "\n\n", sep="")



	##
	##Check arguments, which should be passed in the form <NAME_1>=<VALUE> <NAME_2>=<VALUE> etc
	##
	if( length(my.args) == 0 )
	{
	msg <- paste("This script should be invoked with the following (named) arguments:\n\t",
					paste(EXPECTED.PARAMS, collapse="\n\t"), "\n", sep="")
	stop(msg)
	}
	



cat("Attempting to load required R libraries:\n")

require(QDNAseq)	##this _is_ installed for module load R/3.5.1-foss-2016b-BABS
require(Biobase)
require(BiocGenerics)	##for the magic function "combine"




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
FLAG.THRESHOLD <- MINIMAL.COVERAGE.FLAG.THRESHOLD


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


EXPEC <- list()
	for(i in seq_along(tmp[,1]))
	{
	EXPEC[[ tmp[i,1] ]] <- tmp[i,2]
	}

load(EXPEC[["PARAMS_LOOKUP_FILE"]])
PARAMS[["PARAMS_LOOKUP_FILE"]] <- EXPEC[["PARAMS_LOOKUP_FILE"]]
PARAMS[["CONFIG.R_PATH"]] <- config.path
PARAMS[["COMPLETION_FILE_NAME"]] <- EXPEC[["COMPLETION_FILE_NAME"]]

	##
	##Add other.params if any have been implemented
	##

	###########################################
	###	Compile the alignment metrics	###
	###########################################

cat("\tAttempting to merge alignment metrics:\n")
tmp <- dir(PARAMS[["BAM_SORTED_DIR"]])

regexp <- "_sorted_metrics\\.txt$"

ind <- grep(regexp, tmp)
L <- length(ind)
	if( L == 0 )
	{
	msg <- paste("Failed to detect any matches to \"",regexp , "\" in the sorted bam directory:\n\t\"", PARAMS[["BAM_SORTED_DIR"]], "\"\n", sep="")
	stop(msg)
	}

SUMMARIES <- list()

REPORT <- vector()

	for(i in seq_along(ind) )
	{
	this.ind <- ind[i]
	this.file <- tmp[this.ind]

	full.path <- paste(PARAMS[["BAM_SORTED_DIR"]], this.file, sep="")
	this.sample <- gsub(regexp, "", this.file)


	dat <- as.matrix(read.delim(full.path, sep="\t", header=TRUE))

	SUMMARIES[[ this.sample ]] <- dat


	add <- c(this.sample, dat)
	REPORT <- rbind(REPORT, add)

	}
CN <- c( "Sample_ID", colnames(SUMMARIES[[1]]))
colnames(REPORT) <- CN


	##
	##Replace the full paths to bam files with the raw sample names.
	##

DESIGN <- PARAMS[[ "DESIGN"]]
M <- match(REPORT[, "Sample_ID"], DESIGN[, "LIMS_ID"])
m <- which(is.na(M))
	if( length(m)> 0 )
	{
	msg <- paste("The following LIMS IDs were not found in the look-up table:\n\t",
			paste(REPORT[,"Sample_ID"][m], collapse="\n\t"),
			"\n", sep="")
	stop(msg)

	}

INIT.REPORT <- REPORT
colnames(INIT.REPORT)[colnames(REPORT) == "Sample"] <- "Bam_file"
REPORT <- cbind( DESIGN[, "SAMPLE_NAME_RAW"][M], REPORT[,1], DESIGN[, "GENOME_TAG"][M], REPORT[,-c(1,2), drop=FALSE] )
CN <- colnames(INIT.REPORT)
INIT.REPORT <- cbind(INIT.REPORT, DESIGN[, "ANNO_BED"][M])
colnames(INIT.REPORT) <- c( CN, "ANNO_BED")

#REPORT <- cbind(CLARITY.LOOKUP[, "CLARITY_IDS"][M], REPORT[,1], CLARITY.LOOKUP[, "GENOME_TAGS"][M], REPORT[,-1,drop=FALSE])
colnames(REPORT) <- c("Sample_Name", "LIMS_ID", "Genome_Tag", colnames(SUMMARIES[[1]])[-1])

FLAG <- rep("", dim(REPORT)[1])
aln.per <- as.numeric(gsub(" +", "", REPORT[, "Percent_Aln"]))

w <- which(aln.per < FLAG.THRESHOLD )
	if( length(w) > 0 )
	{
	FLAG[w] <- "<="
	}
CN <- colnames(REPORT)
REPORT <- cbind(REPORT, FLAG)
colnames(REPORT) <- c( CN, "FLAGS")
row.names(REPORT) <- NULL

##
##Add a genome build column following the Genome Tag col
##

genome.tab <- as.matrix(read.delim(GENOMES.LOOKUP.FILE, sep="\t", header=TRUE))

M <- match( REPORT[, "Genome_Tag"], genome.tab[, "Tag"])
m <- which(is.na(M))

GENOME.BUILD <- rep(NA, length(M))
	if( length(m) > 0 )
	{
	msg <- paste("Out of Coffee Error -- unable to re-match the following genome tags:\n\t",
				paste(REPORT[,"Genome_Tag"][m], collapse="\n\t"),
				"\n", sep="")
	cat(msg)

	GENOME.BUILD[-m] <- genome.tab[, "Description"][M[-m]]

	}else
	{
	GENOME.BUILD <- genome.tab[, "Description"][M]
	}

CN <- colnames(REPORT) 
lim <- which(CN == "Genome_Tag")

REPORT <- cbind(REPORT[, 1:lim], GENOME.BUILD, REPORT[,-c(1:lim)])
colnames(REPORT) <- c( CN[1:lim], "Genome_Build", CN[-c(1:lim)])



report.name <- paste(PARAMS[["RESULTS_DIR"]], PARAMS[["PROJECT_NAME"]], "_SequencingMetrics.txt", sep="")
write.table(REPORT, file=report.name, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




	###########################################################
	###	Set up the detailed reports per species.	###
	###########################################################


these.genomes <- unique(REPORT[, "Genome_Tag"])
	for(i in seq_along(these.genomes))
	{
	this.genome <- these.genomes[i]

	this.ind <- which(REPORT[, "Genome_Tag"] == this.genome )
	these.samples <- REPORT[, "Sample_Name"][this.ind]
	these.bams <- INIT.REPORT[, "Bam_file"][this.ind]

	these.detailed.reports <- gsub("\\.bam$", "_detailed_metrics.txt", these.bams)
	VALS <- list()

	CHR <- vector()

		for(j in seq_along(these.samples))
		{
		this.report <- these.detailed.reports[j]

		raw.report <- as.matrix(read.delim(this.report, sep="\t", header=TRUE))

			if( j == 1 )
			{
			CHR <- raw.report[,1]

			}else
			{


			w <- which(CHR != raw.report[,1])
				if( length(w) > 0 )
				{
				msg <- paste("Annotation clash detected while merging detailed reports for \"", this.genome, "\"\n", sep="")
				stop(msg)
				}
			}

		VALS <- cbind( VALS, gsub(" +", "", raw.report[,2]))
		}

	VALS <- cbind(CHR, VALS)
	colnames(VALS) <- c("Chromosome", these.samples)

	detailed.report.name <- paste(PARAMS[["RESULTS_DIR"]], PARAMS[["PROJECT_NAME"]], "_", gsub(" +", "_", this.genome), "_DetailedSequencingMetrics.txt", sep="")
	write.table(VALS, file=detailed.report.name, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

	}



	###################################################
	###	Build the summary QDNAseq figure	###
	###################################################

cat("\tAttempting to merge QDNAseq objects.\n")
tmp <- dir(PARAMS[["QDNASEQ_DIR"]])

tmp.xy <- dir(PARAMS[["XY_QDNASEQ_DIR"]])

regexp <- "_qdnaseq-copy-number\\.RDat$"
ind <- grep(regexp, tmp)

	if( length(ind) == 0 )
	{
	msg <- paste("No matches to \"", regexp, "\" in QDNAseq directory:\n\t\"", PARAMS[["QDNASEQ_DIR"]], "\"\n", sep="")
	stop(msg)
	}

regexp.xy <- "_xy-qdnaseq-copy-number\\.RDat$"
ind.xy <- grep(regexp.xy, tmp.xy)
	if( length(ind.xy) == 0 )
	{
	msg <- paste("No matches to \"", regexp.xy, "\" in XY_QDNAseq directory:\n\t\"", PARAMS[["XY_QDNASEQ_DIR"]], "\"\n", sep="")
	stop(msg) 
	}


	if( length(ind) != length(ind.xy) )
	{

	msg <- paste("Different number of matches found in QDNASeq directory to XY_QDNAseq directory.\n")
		
	stop(msg)
	}

QDAT <- list()
QDAT.XY <- list()


	for(i in seq_along(ind))
	{
	this.ind <- ind[i]
	this.file <- tmp[ this.ind ]
  
	full.file <- paste(PARAMS[["QDNASEQ_DIR"]], this.file, sep="")
	this.sample <- gsub(regexp, "", this.file)

	load(full.file)
    
	QDAT[[ this.sample ]] <- final.copy.numbers
	}

	for(i in seq_along(ind.xy))
	{
	this.ind <- ind.xy[i]
	this.file <- tmp.xy[ this.ind ]
  
	full.file <- paste(PARAMS[["XY_QDNASEQ_DIR"]], this.file, sep="")
	this.sample <- gsub(regexp.xy, "", this.file)

	load(full.file)
    
	QDAT.XY[[ this.sample ]] <- final.copy.numbers
	}


##
##Split reports for different species in the same run
#####################################################

SPECIES.QDATS <- list()
SPECIES.QDATS.XY <- list()
SPECIES.BEDS <- list()

observed.species <- unique(REPORT[, "Genome_Tag"])

	for(i in seq_along(observed.species))
	{
	this.species <- observed.species[i]

	loc.ind <- which(REPORT[, "Genome_Tag"] == this.species)

	SPECIES.QDATS[[ this.species ]] <- list()
	SPECIES.QDATS.XY[[ this.species ]] <- list()
	SPECIES.BEDS[[ this.species ]] <- vector()

		for(j in seq_along(loc.ind) )
		{
		this.sample <- REPORT[, "LIMS_ID"][loc.ind[j]]
		this.bed <- INIT.REPORT[, "ANNO_BED"][loc.ind[j]]
		SPECIES.QDATS[[ this.species ]][[ this.sample ]] <- QDAT[[ this.sample ]]
		SPECIES.QDATS.XY[[ this.species ]][[ this.sample ]] <- QDAT.XY[[ this.sample ]]
		SPECIES.BEDS[[ this.species ]] <- c( SPECIES.BEDS[[ this.species ]], this.bed )
		names(SPECIES.BEDS[[ this.species ]])[length(SPECIES.BEDS[[ this.species ]])] <- this.sample
		}
	}

rm(QDAT)
rm(QDAT.XY)

	for(II in seq_along(SPECIES.QDATS))
	{
	this.species <- names(SPECIES.QDATS)[II]
	QDAT <- SPECIES.QDATS[[ this.species ]]
	QDAT.XY <- SPECIES.QDATS.XY[[ this.species]]
	BED <- SPECIES.BEDS[[this.species]]

		if( length(QDAT) > 1 )
		{

		all.dat <- QDAT[[1]]

			for(j in 2:(length(QDAT)))
			{

			all.dat <- combine( all.dat, QDAT[[j]] )

			}
		dd <- dim(all.dat)

		N <- length(QDAT)	##dd[2]

		}else
		{
		all.dat <- QDAT[[1]]
		N <- 1
		}


	pdf.name <- paste(PARAMS[["RESULTS_DIR"]], PARAMS[["PROJECT_NAME"]],"_", this.species,  "_QDNAseq-copy-number.pdf", sep="")###added this.species here
	pdf.name <- gsub(" ", "_", pdf.name)
	cat("\tPlotting QDNAseq profiles to \"", pdf.name, "\"\n", sep="")
	pdf(pdf.name, width = 2*7, height = (N*7)/2, onefile=TRUE)
	par(mfrow = c(N, 1))

		for(j in seq_along(QDAT) )
		{

		this.bed.file <- BED[j]
			if( !is.na(this.bed.file) )
			{
			cat("Importing annotation bed file:\n\t\"", this.bed.file, "\n", sep="")
			this.bed <- as.matrix(read.delim(this.bed.file, sep="\n", header=FALSE, comment =""))
				if( length(grep("^#", this.bed[1])) == 1 )
				{
				this.bed <- this.bed[-1]
				}
		
			this.bed <- vecsplit(this.bed, split="\t")
			
			}else
			{
			this.bed <- matrix(NA, nrow= 0, ncol = 1)
			}


		plot(QDAT[[ j ]])
		add.bed.anno.to.qdna.plot(qdna.object = QDAT[[j]], bed = this.bed)
		}

		for(j in seq_along(QDAT.XY) )
		{

		this.bed.file <- BED[j]
			if( !is.na(this.bed.file) )
			{
			cat("Importing annotation bed file:\n\t\"", this.bed.file, "\n", sep="")
			this.bed <- as.matrix(read.delim(this.bed.file, sep="\n", header=FALSE, comment =""))
				if( length(grep("^#", this.bed[1])) == 1 )
				{
				this.bed <- this.bed[-1]
				}
		
			this.bed <- vecsplit(this.bed, split="\t")
			
			}else
			{
			this.bed <- matrix(NA, nrow= 0, ncol = 1)
			}

		plot(QDAT.XY[[ j ]])
		add.bed.anno.to.qdna.plot(qdna.object = QDAT[[j]], bed = this.bed)
		}

	dev.off()

	}##added an end to the species loop


	##
	##Call the html report generator
	################################
html.completion.file <- paste(PARAMS[["RESULTS_DIR"]], ".html.done", sep="")
	if( !file.exists(html.completion.file) )
	{

	html.cmd <- paste(BINARY.CALLS[["RSCRIPT"]], " ", WRAPPED.FILES[["CREATE_REPORT.R"]], " ",
		"\"CONFIG.R_PATH=", PARAMS[["CONFIG.R_PATH"]], "\" ",
		"\"PARAMS_LOOKUP_FILE=", PARAMS[[ "PARAMS_LOOKUP_FILE" ]], "\" ",	##we should be able to get everything from here rather than re-constructing everything
		"\"COMPLETION_FILE_NAME=", html.completion.file, "\" ",

		sep="")

	
	S <- system( html.cmd, intern=FALSE)
		if( S != 0 )
		{
		stop("HTML report generation failed.\n")
		}

	}

cat("WRAPPED_do_build_project_summary attempting to mark completion:\n")
cmd <- paste("touch ", PARAMS[["COMPLETION_FILE_NAME"]], sep="")
S <- system(cmd, intern=FALSE)
	if( S != 0 )
	{
	msg <- paste("Ironically, creation of overall successful completion marker failed.\n\t\"", 
		PARAMS[["COMPLETION_FILE_NAME"]], "\"\n", sep="")
	stop(msg)
	}



cat("All done!\n")