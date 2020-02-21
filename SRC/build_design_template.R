

	####################################################################################
	##Script to generate a faux design file indicating the expected format and content##
	####################################################################################



	#########################
	##BUILD.DESIGN.TEMPLATE##
	#########################
########################################################################################################################
build.design.template <- function(
				output.file.name
				)
{


HEADER.CN <- c( "##Description", "##Value" )
HEADER.BLOCK <- rbind(
		##Required entries
		c("##PROJECT_NAME", DEFAULT.PROJECT.NAME),
		c("##FASTQ_DIR",  DEFAULT.FASTQ.SOURCE.DIR),
		c("##OUTPUT_ROOT_DIR",	DEFAULT.OUTPUT.ROOT.DIR),
		##Optional entries

		c("##SCHEDULER_OPTIONS", NA)
		)
L.H <- dim(HEADER.BLOCK)[2]


FOOTER.CN <- c("##LIMS_ID", "##SAMPLE_NAME", "##GENOME_TAG", "##ANNO_BED")
L.F <- length(FOOTER.CN)

FOOTER.RN <- c( "X_1", "X_2", "X_3" )


NC <- max(L.F, L.H)

	if( NC > L.F )
	{
	FOOTER.CN <- c( FOOTER.CN, rep("", NC-L.F) )

	}else if( NC > L.H )
	{

	extra.cols <- NC - L.H
	extra <- matrix("", nrow = dim(HEADER.BLOCK)[1], ncol = extra.cols)
	HEADER.BLOCK <- cbind(HEADER.BLOCK, extra)

	}

FOOTER <- matrix("", nrow = length(FOOTER.RN)+ 1, ncol = NC )


FOOTER[1,] <- FOOTER.CN

FOOTER[2,] <- c("LIMS00001", "My_Sample_1", "hg19", NA)
FOOTER[3,] <- c("LIMS00002", "My_Sample_2", "hg19", NA)
FOOTER[4,] <- c("LIMS00003", "My_Sample_3", "hg19", NA)



OUT <- rbind(HEADER.BLOCK, FOOTER)

write.table(OUT, file=output.file.name, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}##end BUILD.DESIGN.TEMPLATE
########################################################################################################################