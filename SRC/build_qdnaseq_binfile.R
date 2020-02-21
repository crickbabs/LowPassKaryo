library(QDNAseq)
library(Biobase)
library(BSgenome)


	##
	##This will approximate the dataframe required by QDNAseq
	#########################################################


##
##This script will create an RDat object for use with QDNASeq and will require the following
##	A bioconductor bsgenomes object corresponding to your reference sequence (see BSgenome::installed.genomes() & BSgenome::available.genomes() )
##	We don't currently make use of explicit blacklist regions but this could be incorporated in the future.
##	A mappability bed file as produced by, e.g. the encodeproject
##
##mm10_mappability.bed checked on 2019/05/09
##https://www.encodeproject.org/files/ENCFF189JOV/@@download/ENCFF189JOV.bed.gz
##
##GRCh38  checked on 2019/05/09
##https://www.encodeproject.org/files/ENCFF011ZFS/@@download/ENCFF011ZFS.bed.gz
##
##hg19    checked on 2019/05/09
##https://www.encodeproject.org/files/ENCFF496IMC/@@download/ENCFF496IMC.bed.gz
##


	###########################################################
	###							###
	###	RUN SPECIFIC PARAMETERS -- CHANGE THESE!	###
	###							###
	###########################################################

THIS.GENOME <- ##Internal Name for this genome build, e.g. "MM10", in the GENOME.SETS list below
THIS.BIN.SIZE <- 20
THIS.OUTDIR <- ##absolute path to your QDNASeq_bins directory
MEAN.MAPPABILITY <- 97.28	##it would be nice to calculate this specifically



	###########################################################################################################################
	###															###
	###	SET UP SPECIES GENOME BUILD PARAMETERS - TO ADD A NEW SPECIES/GENOME BUILD "JUST" ADD A NEW ENTRY TO THIS LIST	###
	###															###
	###########################################################################################################################

##
##Define the genome build specific settings - these will remain the same regardless of the 
##
GENOME.SETS <- list(

	"HG38" = list(
			bsgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
			mappability.bed = ##absolute path tot he mappabaility bed file for this build,
			chr.to.filter = c("X", "Y")	##QDNASeq recomends the removal of sex chromosomes as they confound the normalisation process
			),


	"MM10" = list(
			bsgenome = "BSgenome.Mmusculus.UCSC.mm10",
			mappability.bed = ##absolute path tot he mappabaility bed file for this build,,
			chr.to.filter = c(	##QDNASeq recomends the removal of sex chromosomes as they confound the normalisation process
						##we also remove floating contigs, etc
						"X", "Y",
						"chrX", "chrY", "chrM",
						"chr1_GL456210_random", "chr1_GL456211_random",
						"chr1_GL456212_random", "chr1_GL456213_random", "chr1_GL456221_random",
						"chr4_GL456216_random", "chr4_GL456350_random", "chr4_JH584292_random",
						"chr4_JH584293_random", "chr4_JH584294_random", "chr4_JH584295_random",
						"chr5_GL456354_random", "chr5_JH584296_random", "chr5_JH584297_random",
						"chr5_JH584298_random", "chr5_JH584299_random", "chr7_GL456219_random",
						"chrX_GL456233_random", "chrY_JH584300_random", "chrY_JH584301_random",
						"chrY_JH584302_random", "chrY_JH584303_random", "chrUn_GL456239",
						"chrUn_GL456359",       "chrUn_GL456360",       "chrUn_GL456366",
						"chrUn_GL456367",       "chrUn_GL456368",       "chrUn_GL456370",
						"chrUn_GL456372",       "chrUn_GL456378",       "chrUn_GL456379",
						"chrUn_GL456381",       "chrUn_GL456382",       "chrUn_GL456383",
						"chrUn_GL456385",       "chrUn_GL456387",       "chrUn_GL456389",
						"chrUn_GL456390",       "chrUn_GL456392",       "chrUn_GL456393",
						"chrUn_GL456394",       "chrUn_GL456396",       "chrUn_JH584304"
					)
			)
	)



                        ###########################################################
                        ###########################################################
                        ###########################################################
                        ###                                                     ###
                        ###                                                     ###
###########################             DO NOT EDIT BELOW THIS LINE             ################################
                        ###                                                     ###
                        ###                                                     ###
                        ###########################################################
                        ###########################################################
                        ###########################################################


	if(! THIS.GENOME %in% names(GENOME.SETS) )
	{

	msg <- paste("Unrecognised genome name \"", THIS.GENOME, "\":\nValid genome names are\n\t",
					paste(names(GENOME.SETS), collapse="\n\t"),
					"\n", sep="")
	stop(msg)
	}

THIS.BSGENOME <- GENOME.SETS[[THIS.GENOME]][[ "bsgenome" ]]
DEFAULT.CHROMOSOMES.TO.FILTER <- GENOME.SETS[[THIS.GENOME]][[ "chr.to.filter" ]]

THIS.MAPPABILITY.BED.FILE <- GENOME.SETS[[ THIS.GENOME ]][[ "mappability.bed" ]]



MEAN.BLACKLIST <- 0.01
MEAN.RESIDUAL <- NA


library(THIS.BSGENOME, character.only = TRUE)

this.bsgenome <- eval(parse(text=THIS.BSGENOME))


OUTPUT.CN <- c(
		"chromosome", "start", "end", "bases", "gc", 
			"mappability", "blacklist", "residual", "use"
		)

	if( !exists("BAK.bins") )
	{
	cat("Creating basic bins:\n")
	flush.console()
	bins <- createBins(bsgenome=this.bsgenome, binSize=THIS.BIN.SIZE)

	BAK.bins <- bins
	}


##so now we have the first section. For now, we fiddle the second section 

	if( !exists( "mapping.bed" ) )
	{
	cat("Importing mappability bed file.\n")
	flush.console()
	mapping.bed <- as.matrix(read.delim(THIS.MAPPABILITY.BED.FILE, sep="\t", header=FALSE))

	mapping.chromosomes <- gsub("^chr", "", mapping.bed[,1])
	mapping.starts <- as.numeric(gsub(" +", "", mapping.bed[,2]))
	mapping.ends <- as.numeric(gsub(" +", "", mapping.bed[,3]))
	}



out.bins <- bins
N <- dim(bins)[1]


mappability <- rep(NA, N)
	for(i in seq_along(mappability))
	{
	this.chromosome <- out.bins$chromosome[i]
	this.start <- out.bins$start[i]
	this.end <- out.bins$end[i]
cat(this.chromosome, ":", this.start, "-", this.end, "\n", sep="")
flush.console()
	ind <- which(
			mapping.chromosomes == this.chromosome &
			! (mapping.ends < this.start) &
			! (mapping.starts > this.end )
			)
			

	L <- length(ind) 

		if( L == 0 )
		{

		mappability[i] <- 0

		}else
		{
		


		interval.starts <- mapping.starts[ind]
		interval.ends <- mapping.ends[ind]

		ord <- order(interval.starts)
		delta <- diff(ord)
		w <- which(delta!=1)
			if(length(w) >0 )
			{
			msg <- paste("The mappability bed file was not ordered.\n",
					"\t(Region: Chromosome ", this.chromosome, ":", this.start, "-", this.end, ")\n", sep="")
					
			stop(msg)
			}
		
		ww <- which(interval.ends[-L] > interval.starts[-1] )
			if( length(ww) > 0 )
			{

			msg <- paste("Overlapping intervals detected in bed file.\n",
					"\t(Region: Chromosome ", this.chromosome, ":", this.start, "-", this.end, ")\n", sep="")
					
			stop(msg)
			}

			if( interval.starts[1] < this.start )
			{
			interval.starts[1] <- this.start
			}

			if( interval.ends[L] > this.end )
			{
			interval.ends[L] <- this.end
			}


		interval.sizes <- interval.ends - interval.starts + 1
		total.cover <- sum(interval.sizes)
		
		full.bin <- this.end - this.start + 1

		per <- 100*total.cover/full.bin
		mappability[i] <- per

			if( per < 0 || per > 100 )
			{

			stop("Mappability calculation error.\n")
			}
	

		}
	cat("\tMappability=", mappability[i], "\n", sep="")

	}

out.bins$"mappability" <- mappability
out.bins$"blacklist" <- BLACK.LIST <- rep(MEAN.BLACKLIST, N)
out.bins$"residuals" <- RESIDUALS <- rep(MEAN.RESIDUAL, N)
out.bins$"use" <- rep(TRUE, N)



annotated.bins <- AnnotatedDataFrame(out.bins,
					varMetadata = data.frame(
						labelDescriptions=c(
							'Chromosome name',
							'Base pair start position',
							'Base pair end position',
							'Percentage of non-N nucleotides (of full bin size)',
							'Percentage of C and G nucleotides (of non-N nucleotides)',
							'Average mappability of 50mers with a maximum of 2 mismatches',
							'Percent overlap with ENCODE blacklisted regions',
							'Median loess residual from 1000 Genomes (50mers)',
							'Whether the bin should be used in subsequent analysis steps'
							),
					row.names=colnames(out.bins)
					)
				)


cat("Bin object now exists.\n")

##write it as an R-object somewhere

output.file.name <- paste(THIS.OUTDIR, "QDNAseq-bin_", THIS.BSGENOME, "_", THIS.BIN.SIZE, "k.RDat", sep="")

default.chromosomes.to.filter <- DEFAULT.CHROMOSOMES.TO.FILTER
this.bin.size <- THIS.BIN.SIZE

save(annotated.bins, default.chromosomes.to.filter, this.bin.size, file=output.file.name)
cat("Bin object written to:\n\t\"", output.file.name, "\"\n", sep="")
