if(!require(optparse)) install.packages("optparse"); library("optparse")
if (!require("ParallelStructure")) install.packages("ParallelStructure", repos="http://R-Forge.R-project.org"); library(ParallelStructure)

option_list = list(
  make_option(c("-j", "--joblist"), type="character", default=NULL, help="joblist file", metavar="character"),
  make_option(c("-c", "--ncpu"), type="integer", default=2, help="Number of CPU cores to be used", metavar="integer"),
  make_option(c("-s", "--structure"), type="character", default=NULL, help="Location of the executable command line STRUCTURE program", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="STRUCTURE input file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="./", help="location of folder to write the output files", metavar="character"),
  make_option(c("-n", "--ninds"), type="integer", default=NULL, help="number of samples in input file", metavar="integer"),
  make_option(c("-l", "--nloci"), type="integer", default=NULL, help="number of loci in input file", metavar="integer"),
  make_option(c("-p", "--ploidy"), type="integer", default=NULL, help="ploidy", metavar="integer")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(any(sapply(opt, is.null))){
	print_help(opt_parser)
	stop("Please provide the requested arguments.")
} 


parallel_structure(
	joblist = opt$joblist,
	n_cpu = opt$ncpu,
	structure_path = opt$structure, 
	infile = opt$input, 
	outpath = opt$output,
	numinds = opt$ninds, 
	numloci = opt$nloci, 
	plot_output = 1, 
	label = 1, 
	popdata = 1, 
	popflag = 1, ## popflag column in the input data
	locdata = 0, 
	phenotypes = 0, 
	markernames = 0, 
	mapdist = 0, 
	onerowperind = 1, ## One row per individual 
	phaseinfo = 0, 
	recessivealleles = 0,  
	phased = 0, 
	extracol = 0, 
	missing = -9, 
	ploidy = opt$ploidy, 
	noadmix = 0, 
	linkage = 0, 
	usepopinfo = 0, 
	locprior = 0,
	inferalpha = 1, 
	alpha = 1, 
	popalphas = 0, 
	unifprioralpha = 1, 
	alphamax = 10, 
	alphapropsd = 0.025, 
	freqscorr = 1,   
	onefst = 0, 
	fpriormean = 0.01, 
	fpriorsd = 0.05, 
	inferlambda = 0, 
	lambda = 1, 
	computeprob = 1,   
	pfromflagonly = 0, 
	ancestdist = 0, 
	startatpopinfo = 0, 
	metrofreq = 10, 
	updatefreq = 1, 
	printqhat = 1, ## Print estimated Q for each individual in separate output file
	revert_convert=0,  
	randomize=1)
