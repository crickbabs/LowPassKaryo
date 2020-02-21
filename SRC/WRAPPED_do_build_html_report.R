

		###################################################################################
		###										###
		###	This script will generate an html file summarising the analysis run 	###
		###		including a description of what the logR plots represent.	###
		###										###
		###################################################################################


my.args <- commandArgs(trailingOnly = TRUE)




EXPECTED.PARAMS <- c(
		"CONFIG.R_PATH",
		"PARAMS_LOOKUP_FILE",
		"COMPLETION_FILE_NAME"
		)
VALID.OTHER.PARAMETERS <- vector()

	if( FALSE )
	{
	##cat("Setting parameters for TESTING.\n")
	##my.args <- c()
	}



d <- date()
this.host <- system("hostname", intern = TRUE)

MSG <- paste("\n\nWRAPPED_do_build_report_html.R invoked on \"", this.host, "\"\n\t@", d, "\n\n", sep="")
MSG <- paste(MSG, "\n\nCalled with the following parameters:\n\t", paste(my.args, collapse="\n\t"), "\n\n", sep="")
cat(MSG)



	#############################################################################################
	##Parse arguments, which should be passed in the form <NAME_1>=<VALUE> <NAME_2>=<VALUE> etc##
	#############################################################################################

	if( length(my.args) == 0 )
	{
	msg <- paste("This script should be invoked with the following (named) arguments:\n\t",
					paste(EXPECTED.PARAMS, collapse="\n\t"), "\n", sep="")
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


	for(i in seq_along(EXPECTED.PARAMS))
	{
	this.param <- EXPECTED.PARAMS[i]
	this.value <- tmp[M[i],2]

	EXPEC[[ this.param ]] <- this.value
	}

load(EXPEC[["PARAMS_LOOKUP_FILE"]])
PARAMS[["COMPLETION_FILE_NAME"]] <- EXPEC[[ "COMPLETION_FILE_NAME" ]]


	###############################
	##Initialise the report file.##
	###############################


report.out <- paste("<HTML>\n<HEAD>\n<TITLE>\n", PARAMS[["PROJECT_NAME"]], " AutoKaryotyping Report\n", "</TITLE>\n</HEAD>\n", sep="")
report.out.file <- paste( PARAMS[["RESULTS_DIR"]],PARAMS[["PROJECT_NAME"]], "_report.html", sep="")



report.out <- paste(report.out, "<BODY>\n", sep="")


report.out <- paste(report.out, "<CENTER><H1>\nAuto-Karyotyping Report for \"", PARAMS[["PROJECT_NAME"]], "\"\n</H1></CENTER>\n", sep="")


report.out <- paste(report.out, 
			"<P><CENTER><H3>", 
			"If you have any queries relating to these data, please contact ", CONTACT.ADDRESS, 
			" who will be able to put you in touch with best person(s) to help you.",
			"</H3></CENTER><P>\n", 
			sep="")

	##
	##We will need
	##

	##sequencing QC report
	##description of what QDNAseq does and how
	##the QDNAseq plots
	##version info


	#############
	##Import QC##
	#############

seq.qc.file <- paste(PARAMS[["RESULTS_DIR"]], PARAMS[["PROJECT_NAME"]], "_SequencingMetrics.txt", sep="")
	if( !file.exists(seq.qc.file) )
	{
	msg <- paste("Could not find sequencing QC file:\n\t\"", seq.qc.file, "\"\n", sep="")
	stop(msg)

	}

qc.table <- as.matrix(read.delim(seq.qc.file, sep="\t", header=TRUE))



report.out <- paste(report.out, "<CENTER>\n", sep="")
report.out <- paste(report.out, "<TABLE BORDER=1>\n\t<TR>\n\t\t", sep="")



	for(i in seq_along(colnames(qc.table)))
	{

	report.out <- paste(report.out, "<TD>", colnames(qc.table)[i], "</TD>", sep="")

	}
report.out <- paste(report.out, "\n\t</TR>\n", sep="")

	for(i in seq_along(qc.table[,1]))
	{
	report.out <- paste(report.out, "\t<TR>\n\t\t", sep="")

		for(j in seq_along(qc.table[i,]))
		{
		report.out <- paste(report.out, "<TD>", qc.table[i,j], "</TD>", sep="")
		}
	report.out <- paste(report.out, "\n\t</TR>\n", sep="")
	}


report.out <- paste(report.out, "</TABLE></CENTER>\n", sep="")





	##
	##Summary of the QDNAseq method here
	##


report.out  <- paste(report.out, "\n<BR>\n<H2>Description</H2>\n", sep="")


report.out <- paste(report.out, "<BR>\n", 
		"Following alignment to the reference genome, copy number estimation was performed using the QDNASeq package:\n",
		"<LIST><LI>",
		"<A HREF=\"https://www.ncbi.nlm.nih.gov/pubmed/25236618\">Scheinin et al \"DNA copy number analysis of fresh and formalin-fixed specimens by shallow whole-genome sequencing with identification and eclusion of problematic regions in the genome assembly.\"</A>\n",
		"</LIST>",
		"<P><BR>In brief, the genome is subdivided into bins of fixed width (1000kb by default) and the number of reads mapping within each bin is calculated.\n",
		"<BR>These \"raw\" counts are then corrected for local GC content and mapability by estimating the median count across bins of the same GC content and mapability.\n",
		"<BR>\"Smoothing\" is then performed by fitting a LOESS surface through the medians.\n",
		"<BR>Finally, the raw counts are corrected by dividing the count for each bin by the LOESS fit corresponding to its GC and mapability.\n",
		"<BR>The plots in the associated PDF files display the log2 transformed ratios.\n",
		"<P><BR>One might think of this process as using each sample to estimate its own \"expected\" number of reads for a given GC-content and mapability conbination, with the smoothing\n",
		"<BR>mitigating outliers based on the assumption that regions of identical GC-content and similar mapability should have similar counts, and simultaneously,\n",
		"<BR>regions of similar GC-content and identical mapability should also have similar counts.\n",
		"<BR>Presentation on the log2 ratios of raw counts to \"expected\" then allows for the standard interpretation: a log2 ratio of 1 indicates that that specific bin has twice as many reads as expected given the other bins for <I>that sample</I>\n",
		"<BR>It should be noted then that these plots do not necessarily indicate \"absolute\" copy number estimates -- for example, in the event that the entire genome is doubled, all regions will be equally over-represented, leading to log2 fold changes of zero rather than 1. Further, the absolute values of the log fold change estimates are <I><B>not</B></I> directly comparable between samples.\n",
		"<P><BR>Finally: the sex chromosomes are explicitly excluded from the median and smoothing process described above.\n",
		##"<BR>The plots below similarly exclude them but separate analysis, re-incorporating them after the calculation of the normalisation factors is also performed by default and these results are available on request.\n",
		"\n\n\n", sep="")




report.out <- paste(report.out, "\n<HR><P>\n", sep="")



if( FALSE )
{##why do we need this?
	##
	##Add in the detailed alignment metrics here, one table per species
	##

these.genomes <- unique(qc.table[, "Genome_Tag"])
	for(i in seq_along(these.genomes))
	{
	this.genome <- these.genomes[i]
	this.ind <- which(qc.table[, "Genome_Tag"] == this.genome )
	this.table <- qc.table[this.ind, , drop=FALSE]

#
#How do we link to the metrics tables?
#
	these.bams <- this.table[, "LIMS_ID"]
	these.detailed.metrics <- gsub("\\.bam$", "_detailed_metrics.txt", these.bams)

	DET <- list()

		for(j in seq_along(these.detailed.metrics) )
		{
		this.file <- these.detailed.metrics[j]

		this.sample <- this.table[, "Sample_Name"][j]

		DET[[ this.sample ]]  <- as.matrix(read.delim(this.file, sep="\t", header=TRUE))
		}


	##
	##Now, add these to a species specific table
	##

	report.out <- paste(report.out, "\n<CENTER><H2>", this.genome, " Detailed Alignment Reports (IDXStats)</H2></CENTER>\n\n", sep="")

	report.out <- paste(report.out, "<P><CENTER>\n", sep="")
	report.out <- paste(report.out, "<TABLE BORDER=1>\n\t<TR>\n\t\t", sep="")
	report.out <- paste(report.out, "<TD>Chromosome</TD> ", sep="")
		for(j in seq_along(DET))
		{
		report.out <- paste(report.out, "<TD>", names(DET)[j], "</TD>", sep="")
		}

	report.out <- paste(report.out, "\n\t</TR>\n", sep="")
	
	THESE.CHROMOSOMES <- DET[[1]][,1]

		for(chr in seq_along(THESE.CHROMOSOMES))
		{
		this.chromosome <- THESE.CHROMOSOMES[chr]

		report.out <- paste(report.out, "\t<TR>\n\t\t<TD>", this.chromosome, "</TD>", sep="")

			for(j in seq_along(DET))
			{
			report.out <- paste(report.out, "<TD>", DET[[j]][ chr ,2], "</TD>", sep="")
			}

		report.out <- paste(report.out, "\n\t</TR>\n", sep="")
		}

	report.out <- paste(report.out, "\n</TABLE>\n\n", sep="")
	report.out <- paste(report.out, "</CENTER>\n", sep="")
	report.out <- paste(report.out, "\n<P>\n<P>\n", sep="")

	}
}



	
report.out <- paste(report.out, "\n<P>\n<HR>\n<P>\n", sep="")
report.out <- paste(report.out, "<H2>Software Version Info:\n<TABLE>\n",sep="")


	for(i in seq_along(BINARY.NAMES))
	{

	this.tag <- names(BINARY.NAMES)[i]

	this.call <- BINARY.CALLS[[ this.tag ]]
	this.name <- BINARY.NAMES[[ this.tag ]]
	this.ver.cmd <- VERSION.COMMANDS[[ this.tag ]]

	this.cmd <- paste(this.call, " ", this.ver.cmd, sep="")
	msg <- paste("\tCalling version for \"", this.name, "\"\n", sep="")
	cat(msg)
	S  <- system(this.cmd, intern=TRUE)

	this.version <- S
	##clean up the trimgalore format...
	this.version <-  gsub("^ +", "", this.version)
	this.version <- paste(this.version[ this.version != "" ], collapse=" ", sep="")

	report.out <- paste(report.out, "\t<TR><TD>", this.name, "</TD><TD>", this.version, "</TD></TR>\n", sep="")
	}


report.out <- paste(report.out, "</TABLE>\n", sep="")


report.out <- paste(report.out, "\n</BODY>\n</HTML>\n", sep="")
write.table(report.out,  file = report.out.file, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)



cmd <- paste("touch ", PARAMS[["COMPLETION_FILE_NAME"]], sep="")
S <- system(cmd, intern=FALSE)
	if( S != 0 )
	{
	msg <- paste("Somewhat ironically, writing of the successful completion marker failed:\n\tCould not \"touch\"\n\t\"", PARAMS[["COMPLETION_FILE_NAME"]], "\"\n", sep="")
	stop(msg)
	}



cat("HTML Report Done!\n")