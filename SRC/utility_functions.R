
	###########################
	###	PATH.CAP	###
	###########################
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



	###########################
	###	VECSPLIT	###
	###########################
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



	###########################
	###	SINGLESPLIT	###
	###########################
########################################################################################################################
singlesplit <- function(
		X,
		split=" "
		)
{

RET <- matrix(NA, nrow = length(X), ncol = 2)

	for(i in seq_along(X))
	{

	tmp <- unlist(strsplit(X[i], split=split))
	RET[i,1] <- tmp[1]
		if( length(tmp) > 1)
		{
		RET[i,2] <- paste(tmp[-1], collapse=split)
		}
	}


return(RET)
}##end SINGLESPLIT
########################################################################################################################


		##################
		##BANNER.MESSAGE##
		##################
########################################################################################################################
banner.message <- function(
			MSG,
			underline.char = "="
			)
{


MSG <- gsub("\\n$", "", MSG)
N <- nchar(MSG)

msg <- paste("\n", MSG, "\n", paste( rep(underline.char, N), collapse=""), "\n\n", sep="")
cat(msg)
flush.console()

}##end BANNER.MESSAGE
########################################################################################################################



	#################
	##USAGE.MESSAGE##
	#################
########################################################################################################################
usage.message <- function()
{

if(DEBUG){banner.message("usage.message called.")}

ret <- ""

	for(i in seq_along(VALID.COMMAND.LINE.OPTIONS))
	{
	this.option <- names(VALID.COMMAND.LINE.OPTIONS)[i]

	these.vals <- VALID.COMMAND.LINE.OPTIONS[[ this.option ]]

	ret <- paste(ret, "\n\t", this.option, " ", sep="")
		if( !is.na(these.vals$"DEFAULT") )
		{
		ret <- paste(ret, "<", paste(these.vals$DEFAULT, collapse=" "), ">", sep="")
		}
	ret <- paste(ret, "\n\t\t", these.vals$"Description", "\n", sep="")
	}

if(DEBUG){banner.message("Returning from usage.message call.\n")}
return(ret)
}##end USAGE.MESSAGE
########################################################################################################################



		###########################
		##PRINT.AVAILABLE.GENOMES##
		###########################
########################################################################################################################
print.available.genomes <- function(
			genomes.file = GENOMES.LOOKUP.FILE
			)
{
if(DEBUG){banner.message("print.available.genomes called.")}

tab <- as.matrix(read.delim(genomes.file, sep="\t", header=TRUE))

ret <- ""

	for(i in seq_along(tab[,1]))
	{
	this.species <- tab[, "Species"][i]
	this.description <- tab[, "Description"][i]
	this.tag <- tab[, "Tag"][i]

	ret <- paste(ret, "\n\t", this.species, "\t", this.description, "\tTag= ", this.tag, "\n", sep="")
	}

if(DEBUG){banner.message("Returning from print.available.genomes call.\n")}
return(ret)
}##end PRINT.AVAILABLE.GENOMES
########################################################################################################################




		######################
		##PARSE.COMMAND.LINE##
		######################
########################################################################################################################
parse.command.line <- function(
				passed.args
				)
{
if(DEBUG){banner.message("parse.command.line called.")}

##
##this should use the passed.args to populate a list of parameters, or just populate with the defaults for any which are missing
##


	##
	##No args
	##########################
	if( length(passed.args) == 0 )
	{##forced to call usage

	MSG <- paste("No parameters passed to the wrapper script:\n\nValid arguments:\n")
	MSG <- paste(MSG, usage.message(), "\n\n", sep="")
	cat(MSG)
	return(NA)
	}

	##
	##Usage call
	###########################
help.ind <- grep("--help", passed.args)
	if( length(help.ind) > 0 )
	{
	MSG <- paste("Usage message requested: Valid arguments are:\n")
	MSG <- paste(MSG, usage.message(), "\n\n", sep="")
	cat(MSG)
	return(NA)
	}


	##
	##Available genomes
	###########################################
genomes.ind <- grep("--available-genomes", passed.args)
	if( length(genomes.ind) > 0 )
	{
	MSG <- paste("List of available genomes requested: Currently valid genomes (and tags) are:\n", sep="")

	MSG <- paste(MSG, print.available.genomes(GENOMES.LOOKUP.TABLE), "\n\n")

	cat(MSG)
	return(NA)
	}


	##
	##Design template
	#################
template.ind <- grep("--template-design", passed.args)
	if( length(template.ind) > 0 )
	{
	#1) look to see if it's followed by a file name 

	candidate.output.file <- NA

		if( template.ind < length(passed.args) )
		{
		candidate.output.file <- passed.args[template.ind+1]

		g <- grep("^-", candidate.output.file)
			if( length(g) > 0 )
			{
			candidate.output.file <- NA
			
			}else
			{
				if( file.exists(candidate.output.file) )
				{
				MSG <- paste("Provided template file already exists -- if you really wish to write a template to:\n\t\"",
								candidate.output.file, "\"\n",
						"please delete the current file.\n\t(This refusal is to avoid accidentally destroying a real ",
							"design file by passing the wrong switch)\n", sep="")
				stop(MSG)
				}
			}
		}

		if( is.na(candidate.output.file) )
		{
		candidate.output.file <- VALID.COMMAND.LINE.OPTIONS[["--template-design"]][["DEFAULT"]]

		msg <- paste("No (valid) name provided for the template design file: using the default of \"", candidate.output.file, "\"\n", sep="")
		cat(msg)
		}

	build.design.template( candidate.output.file)

	MSG <- paste("Design template file written to:\n\t\"", candidate.output.file, "\"\n", sep="")
	cat(MSG)
	return(NA)
	}


design.ind <- grep("^--design$", passed.args)
	if( length(design.ind) == 0 )
	{
	msg <- paste("No recognised arguments detected: please call with \"--help\" for usage details.\n", sep="")
	stop(msg)
	}

	if( design.ind == length(passed.args) || length(grep("^-", passed.args[design.ind+1])) > 0 )
	{
	msg <- paste("Please provide a design file name following the \"--design\" switch\n",
				"\t-- note that design file names may _not_ begin with a \"-\".\n", sep="")
	stop(msg)
	}

DESIGN.FILE.NAME <- passed.args[design.ind+1]

	if( !file.exists(DESIGN.FILE.NAME) )
	{
	msg <- paste("Couldn't find the provided design file:\n\t\"", DESIGN.FILE.NAME, "\"\n", sep="")
	stop(msg)
	}


if(DEBUG){banner.message("Returning from parse.command.line call.\n")}
return(DESIGN.FILE.NAME)
}##end PARSE.COMMAND.LINE
########################################################################################################################




	###########################
	###	SANITISE.NAMES	###
	###########################
########################################################################################################################
sanitise.names <- function(
			raw.names,##vector of raw input names to be sanitised
			valid.chars = "[A-Za-z0-9-]",
			replacement.char = "-"
			)
{



clean.names <- rep(NA, length(raw.names))
altered.names <- rep(NA, length(raw.names))

	##look for naughty names
	for(i in seq_along(raw.names))
	{
	
	this.raw <- raw.names[i]
	this.clean <- this.raw
	this.altered <- NA

	tmp <- unlist(strsplit(this.raw, split=""))
	ind <- grep(valid.chars, tmp, invert=TRUE)
		if( length(ind) > 0 )
		{
		tmp[ind] <- replacement.char
		
		this.clean <- paste(tmp, collapse="")
		this.altered <- this.clean
		}


	clean.names[i] <- this.clean
	altered.names[i] <- this.altered
	}


ii <- which(!is.na(altered.names))
	if( length(ii) > 0 )
	{
	altered.names[ii] <- paste("\"", raw.names[ii], "\" --> \"", altered.names[ii], "\"", sep="")
	}

altered.names <- altered.names[!is.na(altered.names)]


	if( length(altered.names) > 0 )
	{
	msg <- paste("WARNING: ", length(altered.names), " names were sanitised by \"sanitise.names()\"\n", sep="")
	cat(msg)

	}

	if( length(unique(clean.names)) != length(clean.names) )
	{
	msg <- paste("WARNING: Not all names are unique after sanitising.\n")

	cat(msg)
	
	}


ret <- list(
	"input" = raw.names,
	"clean" = clean.names,
	"altered" = altered.names
	)

return(ret)
}##end SANITISE.NAMES
########################################################################################################################



		##################################
		##	PARSE.DESIGN.FILE	##
		##################################
########################################################################################################################
parse.design.file <- function(
			DESIGN.FILE.NAME,
			GENOME.LOOKUP = GENOME.LOOKUP.TABLE
			)
{
if(DEBUG){banner.message("Parse.design.file called.\n")}

	if( !file.exists(DESIGN.FILE.NAME) )
	{
	msg <- paste("Unable to find provided design file:\n\t\"", DESIGN.FILE.NAME, "\"\n", sep="")
	stop(msg)
	}


raw.design <- as.matrix(read.delim(DESIGN.FILE.NAME, sep="\t", header=FALSE, comment = "", quote=""))



footer.start.row <- grep("^##LIMS_ID", raw.design[,1])
	if( length(footer.start.row) != 1 )
	{
	msg <- paste("Unable lo locate a unique match to \"##LIMS_ID\" in the design file\n\t(\"", DESIGN.FILE.NAME, "\")\n", sep="")
	stop(msg)
	}

	if( footer.start.row == dim(raw.design)[1] )
	{
	msg <- paste("No sample info detected in design file\n\t(\"", DESIGN.FILE.NAME, "\")\n", sep="")
	stop(msg)
	}

	if( footer.start.row == 1 )
	{
	msg <- paste("No header info detected in design file\n\t(\"", DESIGN.FILE.NAME, "\")\n", sep="")
	stop(msg)
	}


FOOTER.IND <- (footer.start.row+1):(dim(raw.design)[1])

FOOTER.BLOCK <- raw.design[ FOOTER.IND,,drop=FALSE]
colnames(FOOTER.BLOCK) <- raw.design[ footer.start.row,]

HEADER.BLOCK <- raw.design[ 1:(footer.start.row-1),,drop=FALSE]
colnames(HEADER.BLOCK) <- NULL

PARAMS <- list()
err.msg <- ""

		###################################
		###				###
		###	Parse the Header block	###
		###				###
		###################################


##PROJECT_NAME##
################
w <- which(HEADER.BLOCK[,1] == "##PROJECT_NAME")
	if( length(w) != 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Failed to detect a unique match to \"##PROJECT_NAME\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")

	}
PARAMS[["PROJECT_NAME"]] <- HEADER.BLOCK[w,2]


##FASTQ_DIR##
#############
w <- which(HEADER.BLOCK[,1] == "##FASTQ_DIR")
	if( length(w) != 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Failed to detect a unique match to \"##FASTQ_DIR\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")
	}
PARAMS[["FASTQ_DIR"]] <- path.cap(HEADER.BLOCK[w,2])
	if( !file.exists(PARAMS[["FASTQ_DIR"]]) )
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate provided FastQ source directory:\n\t\"##FASTQ_DIR\" = \"", PARAMS[["FASTQ_DIR"]], "\"\n\n", sep="")

	}


##OUTPUT_ROOT_DIR##
###################
w <- which( HEADER.BLOCK[,1] == "##OUTPUT_ROOT_DIR" )
	if( length(w) != 1)
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to \"##OUTPUT_ROOT_DIR\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")
	}
PARAMS[["OUTPUT_ROOT_DIR"]] <- path.cap( HEADER.BLOCK[w,2] )
	if( !file.exists(PARAMS[["OUTPUT_ROOT_DIR"]] ) )
	{
	
	cmd <- paste("mkdir ", PARAMS[[ "OUTPUT_ROOT_DIR" ]], sep="")
	S <- system(cmd, intern=FALSE)
		if( S != 0 )
		{
		msg <- paste("Attempt to create output root directory failed:\n\t\"##OUTPUT_ROOT_DIR\" = \"", PARAMS[["OUTPUT_ROOT_DIR"]], "\"\n\n", sep="")
		stop()
		}
	}

##SCHEDULER_OPTIONS##
#####################

w <- which( HEADER.BLOCK[,1] == "##SCHEDULER_OPTIONS" )
	if( length(w) > 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to \"##SCHEDULER_OPTIONS\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")

	}else if( length(w) == 0 || is.na(HEADER.BLOCK[w,2]))
	{
	PARAMS[["SCHEDULER_OPTIONS"]] <- DEFAULT.SCHEDULER.OPTIONS

	}else
	{
	PARAMS[["SCHEDULER_OPTIONS"]] <- HEADER.BLOCK[w,2]

	}


		###################################
		###				###
		###	Parse the Footer block	###
		###				###
		###################################


F.CN <- colnames(FOOTER.BLOCK)

LIMS.COL <- which(F.CN == "##LIMS_ID")
SAMPLE.COL <- which(F.CN == "##SAMPLE_NAME" )
GENOME.COL <- which( F.CN == "##GENOME_TAG" )
ANNO.COL <- which(F.CN == "##ANNO_BED" )

	if( length(LIMS.COL) != 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to \"##LIMS_ID\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")
	}
	if( length(SAMPLE.COL) != 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to \"##SAMPLE_NAME\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")
	}
	if( length(GENOME.COL) != 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to \"##GENOME_TAG\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")
	}
	if( length(ANNO.COL) > 1 )
	{
	err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to \"##ANNO_BED\" in design file.\n\t(\"", DESIGN.FILE.NAME, "\")\n\n", sep="")
	}

LIMS.IDS <- FOOTER.BLOCK[, LIMS.COL]
SAMPLE.NAME.RAW <- FOOTER.BLOCK[, SAMPLE.COL]
GENOME.TAGS <- FOOTER.BLOCK[, GENOME.COL]

	if( length(ANNO.COL) == 0 )
	{
	ANNO.BEDS <- rep(DEFAULT.ANNO.BED.FILE, length(LIMS.IDS))

	}else
	{
	ANNO.BEDS <- FOOTER.BLOCK[, ANNO.COL]
	}


	##
	##Check we can find matching FastQ files
	########################################
d <- dir( PARAMS[["FASTQ_DIR"]] )
	for(i in seq_along(LIMS.IDS))
	{
	this.lims.id <- LIMS.IDS[i]

	g <- grep(paste("^", this.lims.id, "_", sep=""), d)
		if( length(g) == 0 )
		{
		err.msg <- paste(err.msg, ">>>>>Unable to locate a  match to LIMS_ID \"",this.lims.id ,"\" in the provided FastQ directory.\n\t(\"", PARAMS[["FASTQ_DIR"]], "\")\n\n", sep="")
		}
	}


	##
	##Check we can find the genome tags
	###################################

	for( i in seq_along(GENOME.TAGS))
	{
	this.tag <- GENOME.TAGS[i]
	w <- which(GENOME.LOOKUP[, "Tag"] == this.tag)
		if( length(w) != 1 )
		{
		err.msg <- paste(err.msg, ">>>>>Unable to locate a unique match to the genome tag \"", this.tag, "\" associated with LIMS_ID \"", LIMS.IDS[i], "\" in the design file.\n", sep="")
		}
	}

	##
	##Finaly, check we can find any annotation bed files which may have been provided.
	##################################################################################


	for(i in seq_along(ANNO.BEDS))
	{
	this.bed <- ANNO.BEDS[i]

		if( !is.na(this.bed) && !file.exists( this.bed ))
		{
		err.msg <- paste(">>>>>Unable to locate annotation bed file associated with the LIMS_ID \"", LIMS.IDS[i], "\":\n\t\"", this.bed, "\"\n", sep="")
		}
	}



DESIGN.CN <- c( "LIMS_ID", "SAMPLE_NAME_RAW", "SAMPLE_NAME_SANITISED", "GENOME_TAG", "ANNO_BED")



TMP <-  sanitise.names(
			raw.names = SAMPLE.NAME.RAW
			)

#########
#sanitise.names
#########
#ret <- list(
#	"input" = raw.names,
#	"clean" = clean.names,
#	"altered" = altered.names
#	)
#########

SAMPLE.NAME.CLEAN <- TMP$"clean"
DESIGN <- cbind(
		LIMS.IDS, 
		SAMPLE.NAME.RAW, 
		SAMPLE.NAME.CLEAN, 
		GENOME.TAGS, 
		ANNO.BEDS
		)
colnames(DESIGN) <- DESIGN.CN
PARAMS[["DESIGN"]] <- DESIGN

		########
		##DONE##
		########

	if( nchar(err.msg ) > 0 )
	{
	cat("\n\n\n")
	stop(err.msg)
	}


if(DEBUG){banner.message("Returning from parse.design.file call.\n")}
return(PARAMS)
}##end PARSE.DESIGN.FILE
########################################################################################################################



########################################################################################################################
associate.fastq.with.lims.id <- function(
					params = PARAMS
					)
{
if(DEBUG){banner.message("associate.fastq.with.lims.id called.\n")}

this.design <- params[["DESIGN"]]
LIMS.IDS <- this.design[, "LIMS_ID"]
this.fastq.dir <- params[["FASTQ_DIR"]]



d <- dir(this.fastq.dir)
candidate.lims.ids <- vecsplit(d, split="_")

FASTQ.FILES <- matrix(NA, nrow = 0, ncol = 2)

lims.id.errs <- ""


r2.ind <- grep(LOCAL.R2.REGEXP,	d)
	if( length(r2.ind) == 0 )
	{
	cat("This dataset is inferred to be Single End only.\n")
	SINGLE.ENDED <- TRUE

	}else
	{
	cat("This dataset is inferred to be Paired End only.\n")
	SINGLE.ENDED <- FALSE

	}

	for(i in seq_along(LIMS.IDS))
	{

	this.lims.id <- LIMS.IDS[i]
	this.regexp <- paste("^", this.lims.id, "_", sep="")
	g <- grep(this.regexp, d)

	L <- length(g)
		##
		##
		##

		if( L == 0 )
		{
		msg <- paste("No matches to the LIMS ID \"", this.lims.id, "\" in the provided fastq directory.\n", sep="")
		#cat(msg)

		lims.id.errs <- paste( lims.id.errs, msg, sep="")
		next()
		}

		if( L == 1 && !SINGLE.ENDED  )
		{
		msg <- paste("Single matche to the LIMS ID \"", this.lims.id, "\" in the provided fastq directory but run inferred to be paried end\n", sep="")	
		lims.id.errs <- paste( lims.id.errs, msg, sep="")
		next()
		}

		if( L == 2 && SINGLE.ENDED  )
		{
		msg <- paste("Two matches to the LIMS ID \"", this.lims.id, "\" in the provided fastq directory but run inferred to be single ended\n", sep="")	
		lims.id.errs <- paste( lims.id.errs, msg, sep="")
		next()
		}

		if( L > 2 )
		{
		msg <- paste("Multiple (", L, ") matches to LIMD ID \"", this.lims.id, "\" in the provided fastq directory.\n", sep="")
		lims.id.errs <- paste( lims.id.errs, msg, sep="")
		next()
		}


	these.files <- d[g]
	r1.ind <- grep(LOCAL.R1.REGEXP, these.files)##we're explicitly assuming the illumina pipeline naming remains fixed
		if( length(r1.ind) != 1 )
		{
		msg <- paste("Failed to find a unique R1 match for LIMS ID \"", this.lims.id, "\" (found ", 
		length(r1.ind), ") in the provided fastq directory.\n", sep="")

		cat(msg)
		lims.id.errs <- paste( lims.id.errs, msg, sep="")
		next()
		}



		if( SINGLE.ENDED )
		{
		FASTQ.FILES <- rbind( FASTQ.FILES, 
					c( these.files[ r1.ind ], "SINGLE_END" )
					)

		}else##paired ended!
		{

			if( L == 1 )
			{
			msg <- paste("Only one match to the LIMS ID \"", this.lims.id, 
					"\" in the provided fastq directory (paired end data expected)\n", sep="")
			cat(msg)
			lims.id.errs <- paste( lims.id.errs, msg, sep="")
			next()
			}

			if( L > 2 )	
			{	
			msg <- paste("More matches to LIMS ID \"", this.lims.id, "\" found (", L, 
					") than expected (",2,") in the provided fastq directory.\n", sep="")

			cat(msg)
			lims.id.errs <- paste( lims.id.errs, msg, sep="")
			}

		r2.ind <- grep(LOCAL.R2.REGEXP, these.files)	

			if( length(r2.ind) != 1 )
			{
			msg <- paste("Failed to find a unique R2 match for LIMS ID \"", this.lims.id, "\" (found ", 
						length(r2.ind), ") in the provided fastq directory.\n", sep="")
	
			cat(msg)
			lims.id.errs <- paste( lims.id.errs, msg, sep="")
			next()
			}


		FASTQ.FILES <- rbind( FASTQ.FILES, 
				c( these.files[ r1.ind ], these.files[ r2.ind ] )
				)

		}

	}

	if( nchar(lims.id.errs) > 0 )
	{
	stop(lims.id.errs)
	}


return(FASTQ.FILES)

if(DEBUG){banner.message("Returning from associate.fastq.with.lims.id call.\n")}
}##end ASSOCIATE.FASTQ.WITH.LIMS.ID
########################################################################################################################





	#############################
	##ADD.BED.ANNO.TO.QDNA.PLOT##
	#############################
########################################################################################################################
add.bed.anno.to.qdna.plot <- function(
				qdna.object,
				bed
				)
{


##
##Identify the chromosome offsets
##

RN <- row.names(fData(qdna.object))
TMP <- vecsplit(RN, split=":")
these.chrs <- TMP[,1]
TMP.2 <- vecsplit(TMP[,2], split="-")
these.starts <- as.numeric(TMP.2[,1])
these.ends <- as.numeric(TMP.2[,2])

obs.chrs <- unique(these.chrs)

w <- which(!bed[,1] %in% obs.chrs)
	if( length(w) > 0 )
	{
	msg <- paste("The following chromosomes from the bed annotation file were not recognised by the QDNASeq object:\n\t",
						paste(bed[w,1], collapse="\n\t"),
						"\n", 
						sep="")
	stop(msg)
	}	


start.indices <- rep(NA, length(obs.chrs))
	for(i in seq_along(obs.chrs))
	{
	this.chr <- obs.chrs[i]
	ind <- which( these.chrs == this.chr )

	start.indices[i] <- min(ind)

	}

end.indices <- c( 1, start.indices[-1]-1)



START.POINTS <- cumsum(these.ends[ end.indices] )
names(START.POINTS) <- obs.chrs


M <- max(qdna.object@assayData$"copynumber", na.rm=TRUE)
M <- M-1

	if( M < 0 )
	{
	cat("Unable to infer appropriate y-range for annotation.\n")
	return()
	}

	for(i in seq_along(bed[,1]))
	{
	this.entry <- this.bed[i,]

	this.chr <- as.character(this.entry[1])
	this.start <- as.numeric(this.entry[2])
	this.end <- as.numeric(this.entry[3])
	this.anno <- this.entry[4]

	this.offset <- as.numeric(START.POINTS[ as.character(this.chr) ])

	this.start <- this.start + this.offset
	this.end <- this.end + this.offset



	lines(x = rep(this.start, 2), y = c(-M,M), col=DEFAULT.ANNO.COLOUR, lwd=2)
	lines(x = rep(this.end,2), y = c(-M,M), col=DEFAULT.ANNO.COLOUR, lwd=2)

	text( x=(this.end+this.start)/2, y=1, lab=this.anno, pos=3)
	}
}##end ADD.BED.ANNO.TO.QDNA.PLOT
########################################################################################################################


