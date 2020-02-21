


                ###########################
                ###                     ###
                ###     DO_TRIMGALORE   ###
                ###                     ###
                ###########################

        ###########################################################################
        ##This script is intended to be automatically invoked by a wrapper script##
        ###########################################################################



my.args <- commandArgs( trailingOnly = TRUE )

EXPECTED.PARAMS <- c(
                "QCTRIM",           ##system specific method for invoking TRIMGALORE (e.g. including module load statements on camp)
		"INPUT.FASTQ.FILE_R1",	##Use this option to pass single end data
		"INPUT.FASTQ.FILE_R2",	## _IF_ data is single ended, pass NA here _AND_ set the PAIRED.END param to FALSE
		"PAIRED.END",		##Set TRUE if data is paired end and FALSE if not
		"OUTPUT.DIRECTORY",
		"COMPLETION.FILE.NAME"
		)



VALID.OTHER.PARAMS <- list(
	"--quality" = 20,
        		##Trim low-quality ends from reads in addition to adapter removal. 
			##For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. 
			##Other files are quality and adapter trimmed in a single pass. 
			##The algorithm is the same as the one used by BWA (Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; 	
			##	cut sequence at the index at which the sum is minimal).
			##Default Phred score: 20
	"--phred33" = TRUE,
		##Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding) for quality trimming.
		##Default: ON
	"--phred64" = FALSE,
		##Instructs Cutadapt to use ASCII+64 quality scores as Phred scores (Illumina 1.5 encoding) for quality trimming.
	"--fastqc" = TRUE,
		##Run FastQC in the default mode on the FastQ file once trimming is complete.
	"--fastqc_args" = NULL,## "<ARGS>"
		##Passes extra arguments to FastQC. If more than one argument is to be passed to FastQC they must be in the form arg1 arg2 [..].
		##An example would be: --fastqc_args "--nogroup --outdir /home/".
		##Passing extra arguments will automatically invoke FastQC, so --fastqc does not have to be specified separately.
	"--adapter" = NULL, ## <STRING>
        	##Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will try to auto-detect whether the Illumina universal, 
		##Nextera transposase or Illumina small RNA adapter sequence was used. Also see --illumina, --nextera and --small_rna.
        	##If no adapter can be detected within the first 1 million sequences of the first file specified Trim Galore defaults to --illumina. 
		##A single base may also be given as e.g. -a A{10}, to be expanded to -a AAAAAAAAAA.

	"--adapter2" = NULL,	## <STRING>
        	##Optional adapter sequence to be trimmed off read 2 of paired-end files. 
		##This option requires --paired to be specified as well. 
		##If the libraries to be trimmed are smallRNA then a2 will be set to the Illumina small RNA 5' adapter automatically (GATCGTCGGACT). 
		##A single base may also be given as e.g. -a2 A{10}, to be expanded to -a2 AAAAAAAAAA.
	"--illumina" = FALSE,
		##Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter AGATCGGAAGAGC instead of the default auto-detection of adapter sequence.
	"--nextera" = FALSE,
		##Adapter sequence to be trimmed is the first 12bp of the Nextera adapter CTGTCTCTTATA instead of the default auto-detection of adapter sequence.
	"--small_rna" = FALSE,
		##Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter TGGAATTCTCGG instead of the default auto-detection of adapter sequence.
		##Selecting to trim smallRNA adapters will also lower the --length value to 18bp. 
		##If the smallRNA libraries are paired-end then -a2 will be set to the Illumina small RNA 5' adapter 
		##automatically (GATCGTCGGACT) unless -a 2 had been defined explicitly.
	"--max_length" = NULL,	## <INT>
		##Discard reads that are longer than bp after trimming. This is only advised for smallRNA sequencing to remove non-small RNA sequences.
	"--stringency" = 1,	# <INT>
		##Overlap with adapter sequence required to trim a sequence.
		##Defaults to a very stringent setting of 1, i.e. even a single base pair of overlapping sequence will be trimmed of the 3' end of any read.
	"-e" = 0.1,## <ERROR RATE>
		##Maximum allowed error rate (no. of errors divided by the length of the matching region)
		##Default: 0.1
	"--gzip" = FALSE,
		##Compress the output file with gzip.
		##If the input files are gzip-compressed the output files will be automatically gzip compressed as well.
	"--dont_gzip" = FALSE,
		##Output files won't be compressed with gzip. This overrides --gzip.
	"--length" = 20,## <INT>
		##Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of 0 effectively disables this behaviour.
		##Default: 20 bp.
		##For paired-end files, both reads of a read-pair need to be longer than bp to be printed out to validated paired-end files 
		##(see option --paired). If only one read became too short there is the possibility of keeping such unpaired single-end reads (see --retain_unpaired).
		##Default pair-cutoff: 20 bp.
	"--max_n" = NULL,## COUNT
		##The total number of Ns (as integer) a read may contain before it will be removed altogether.
		##In a paired-end setting, either read exceeding this limit will result in the entire pair being removed from the trimmed output files.
	"--trim-n" = FALSE,
		##Removes Ns from either side of the read.
		##This option does currently not work in RRBS mode.
	"--output_dir" = NULL,	##<DIR>
		##If specified all output will be written to this directory instead of the current directory. If the directory doesn't exist it will be created for you.
	"--no_report_file" = FALSE,
		##If specified no report file will be generated.
	"--suppress_warn" = FALSE,
		##If specified any output to STDOUT or STDERR will be suppressed.
	"--clip_R1" = NULL,## <int>
		##Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads). 
		##This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end.
		##Default: OFF
	"--clip_R2" = NULL, ## <int>
		##Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only). 
		##This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end.
		##For paired-end BS-Seq, it is recommended to remove the first few bp because the end-repair reaction may 
		##introduce a bias towards low methylation. Please refer to the M-bias plot section in the Bismark User Guide for some examples.
		##Default: OFF
	"--three_prime_clip_R1" = NULL, #<int>
		##Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end reads) AFTER adapter/quality trimming has been performed. 
		##This may remove some unwanted bias from the 3' end that is not directly related to adapter sequence or basecall quality.
		##Default: OFF
	"--three_prime_clip_R2" = NULL,## <int>
		##Instructs Trim Galore to re move <int> bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed. 
		##This may remove some unwanted bias from the 3' end that is not directly related to adapter sequence or basecall quality.
        	##Default: OFF
	"--2colour" = NULL,
	##A.K.A
	"--nextseq" = NULL, 
        	##This enables the option --nextseq-trim=3'CUTOFF within Cutadapt, 
		##which will set a quality cutoff (that is normally given with -q instead), 
		##but qualities of G bases are ignored. This trimming is in common for the NextSeq- and NovaSeq-platforms, 
		##where basecalls without any signal are called as high-quality G bases. 
		##More on the issue of G-overcalling may be found here: 
		##https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/. 
		##This is mutually exlusive with -q INT.
	"--path_to_cutadapt" = NULL,## </path/to/cutadapt>
		##You may use this option to specify a path to the Cutadapt executable, e.g. /my/home/cutadapt-1.7.1/bin/cutadapt. 
		##Else it is assumed that Cutadapt is in the PATH.
	"--basename" = NULL,## <PREFERRED_NAME>
		##Use PREFERRED_NAME as the basename for output files, instead of deriving the filenames from the input files. 
		##Single-end data would be called PREFERRED_NAME_trimmed.fq(.gz), 
		##or PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz) for paired-end data. 
		##--basename only works when 1 file (single-end) or 2 files (paired-end) are specified, but not for longer lists.
	"--cores" = 1,## INT
		##Number of cores to be used for trimming [default: 1]. 
		##For Cutadapt to work with multiple cores, it requires Python 3 as well as parallel gzip (pigz) installed on the system. 
		##The version of Python used is detected from the shebang line of the Cutadapt executable (either cutadapt, or a specified path). 
		##If Python 2 is detected, --cores is set to 1. If pigz cannot be detected on your system, Trim Galore reverts to using gzip compression. 
		##Please note that gzip compression will slow down multi-core processes so much that it is hardly worthwhile, 
		##please see: https://github.com/FelixKrueger/TrimGalore/issues/16#issuecomment-458557103 for more info).
		##########
		##Actual core usage: It should be mentioned that the actual number of cores used is a little convoluted. 
		##Assuming that Python 3 is used and pigz is installed, --cores 2 would use 2 cores to read the input (probably not at a high usage though), 
		##2 cores to write to the output (at moderately high usage), and 2 cores for Cutadapt itself + 2 additional cores for Cutadapt 
		##(not sure what they are used for) + 1 core for Trim Galore itself. 
		##So this can be up to 9 cores, even though most of them won't be used at 100% for most of the time. 
		##Paired-end processing uses twice as many cores for the validation (= writing out) step. 
		##--cores 4 would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15, and so forth.
		##
		##It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
			###############################################################
		###########	SPECIFIC TRIMMING - without adapter/quality trimming###########
			###############################################################
	"--hardtrim5" = NULL, ##<int>
		##Instead of performing adapter-/quality trimming, this option will simply hard-trim sequences to bp from the 3'-end. 
		##Once hard-trimming of files is complete, Trim Galore will exit. 
		##Hard-trimmed output files will end in .<int>bp_5prime.fq(.gz).
	"--hardtrim3" = NULL, ##<int>
		##Instead of performing adapter-/quality trimming, this option will simply hard-trim sequences to bp from the 5'-end. 
		##Once hard-trimming of files is complete, Trim Galore will exit. Hard-trimmed output files will end in .<int>bp_3prime.fq(.gz).
	"--clock" = FALSE,
		##See documentation for full details
			###########################################################
		###########	RRBS-specific options (MspI digested material):	###########
			###########################################################
	"--rrbs" = FALSE,
		##Specifies that the input file was an MspI digested RRBS sample (recognition site: CCGG). 
		##Sequences which were adapter-trimmed will have a further 2 bp removed from their 3' end. 
		##This is to avoid that the filled-in C close to the second MspI site in a sequence is used for methylation calls. 
		##Sequences which were merely trimmed because of poor quality will not be shortened further.
	"--non_directional" = FALSE,
		##Selecting this option for non-directional RRBS libraries will screen quality-trimmed sequences for CAA or CGA at the start of the read and, 
		##if found, removes the first two base pairs. 
		##Like with the option --rrbs this avoids using cytosine positions that were filled-in during the end-repair step. 
		##--non_directional requires --rrbs to be specified as well.
	"--keep" = FALSE,
		##Keep the quality trimmed intermediate file. If not specified the temporary file will be deleted after adapter trimming. 
		##Only has an effect for RRBS samples since other FastQ files are not trimmed for poor qualities separately.
		##Default: OFF
		##########	
		##Note for RRBS using MseI:
		##
		##If your DNA material was digested with MseI (recognition motif: TTAA) instead of MspI it is NOT necessary to specify 
		##--rrbs or --non_directional since virtually all reads should start with the sequence TAA, 
		##and this holds true for both directional and non-directional libraries. 
		##As the end-repair of TAA restricted sites does not involve any cytosines it does not need to be treated especially. 
		##Instead, simply run Trim Galore! in the standard, i.e. non-RRBS, mode.
		##
			###########################################
		###########	Paired-end specific options:	###########
			###########################################
	"--paired" = FALSE,
		##This option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files. 
		##To pass the validation test, both sequences of a sequence pair are required to have a certain minimum length which is 
		##governed by the option --length (see above). 
		##If only one read passes this length threshold the other read can be rescued (see option --retain_unpaired).
		##Using this option lets you discard too short read pairs without disturbing the sequence-by-sequence order of 
		##FastQ files which is required by many aligners. 
		##Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g. file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... .
	"--trim1" = FALSE,
		##Trims 1 bp off every read from it's 3' end
		##This may be needed for FastQ file which are to be aligned as paried-end data with bowtie1
		##
	"--retain_unpaired" = FALSE,
		##If only one of the two paired-end reads became too short, the longer read will be written to either 
		##.unpaired_1.fq or .unpaired_2.fq output files. 
		##The length cutoff for unpaired single-end reads is governed by the parameters -r1/--length_1 and -r2/--length_2.
		##Default: OFF.
	"--length_1" = 35,	## <INT>
		##Unpaired single-end read length cutoff needed for read 1 to be written to .unpaired_1.fq output file. 
		##These reads may be mapped in single-end mode.
		##Default: 35 bp
	"--length_2" = 35	## <INT>
		##Unpaired single-end read length cutoff needed for read 2 to be written to .unpaired_2.fq output file. 
		##These reads may be mapped in single-end mode.
		##Default: 35 bp  

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

MSG <- paste("\n\ndo_trimgalore.R invoked on \"", this.host, "\"\n\t@", d, "\n\n", sep="")
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


	if( !file.exists( PARAMS[["INPUT.FASTQ.FILE_R1"]] ) )
	{
	msg <- paste("Unable to locate provided (R1) fastq file:\n\t\"", PARAMS[["INPUT.FASTQ.FILE_R1"]], "\"\n", sep="")
	stop(msg)
	}

	if( PARAMS[["PAIRED.END"]] == "TRUE" )
	{
		if( !file.exists( PARAMS[["INPUT.FASTQ.FILE_R2"]] ) )
		{
		msg <- paste("Unable to locate provided (R2) fastq file:\n\t\"", PARAMS[["INPUT.FASTQ.FILE_R2"]], "\"\n", sep="")
		stop(msg)
		}		
	}


	if( !file.exists( PARAMS[[ "OUTPUT.DIRECTORY" ]]  ) )
	{##all directories to be created _prior_ to calling to avoid sitting in the queue only to fail because the path's wrong.

	msg <- paste("Provided output directory does not exist:\n\t\"", PARAMS[[ "OUTPUT.DIRECTORY"]], "\"\n", sep="")
	stop(msg)
	}



	##############################
	###Build the system command###
	##############################

##USAGE:  'trim_galore [options] <filename(s)>'    or    'trim_galore --help'    for more options


trimgalore.cmd <- paste(
				PARAMS[[ "QCTRIM" ]], " ", 
				" --output_dir ", PARAMS[["OUTPUT.DIRECTORY"]], " ",
				sep="")

	if( PARAMS[["PAIRED.END"]] == "TRUE" )
	{
	trimgalore.cmd <- paste(trimgalore.cmd, " --paired ", sep="")

	}

	for(i in seq_along(other.params[,1]))
	{

	trimgalore.cmd <- paste(trimgalore.cmd, " ", other.params[i,1], " ", other.params[i,2], " ", sep="")
	msg <- paste("Added optional argument \"", other.params[i,1], "\" = \"",other.params[i,2], "\"\n", sep="")
	cat(msg)
	}


trimgalore.cmd <- paste(trimgalore.cmd, " ", PARAMS[["INPUT.FASTQ.FILE_R1"]], " ", sep="")

	if( PARAMS[["PAIRED.END"]] == "TRUE" )
	{
	
	trimgalore.cmd <- paste(trimgalore.cmd, " ", PARAMS[["INPUT.FASTQ.FILE_R2"]], " ", sep="")
	}

cat("Trimgalore command:\n\t\"", trimgalore.cmd, "\"\n", sep="")

S <- system(trimgalore.cmd, intern = FALSE)
	if( S != 0 )
	{
	msg <- "Trimgalore call failed.\n"
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

