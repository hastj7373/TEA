###############################################################
# Functions
###############################################################
# get.calls <- function(ctype, rasym, type, samples=NULL)
###############################################################

source("/bitools/TEA/Tea-0.6.2/R/rid.r")

library(RColorBrewer)

can <- read.table("cancer", header=F)
mat <- read.table("matched", header=F)

nonmatched.controls.hg18 = c()
if(file.exists("nonmatched")){
	nor <- read.table("nonmatched", header=F)
	nonmatched.controls.hg19 = as.vector(nor[,1])
} else{
	nonmatched.controls.hg19 = c()
}

two.col = c("#2B6CA9", "#F25D44") # blue and red
three.col = c("#2B6CA9", "#F25D44", "#3CB371") # blue, red, green

get.calls <- function(ctype, rasym, type, samples=NULL)
{
	pdir = getwd()
	cmd = sprintf("dirname %s", pdir)
	bdir <- system(cmd, intern = TRUE)
	dir = sprintf("%s", pdir)
	
	samples = list.files(dir, "*_(cancer|normal)$")
	names(samples) = samples
	normals = samples[grepl("normal", samples)]
	#print(paste(normals))
	if (is.null(samples)) {
		if (ctype == "brca_trip_negative") {
			samples = c("normal", "primary", "xenograft_from_primary", "brain_metastasis")
			names(samples) = samples
			cancers = tail(samples, -1)
		} else {
			samples = list.files(dir, "*_(cancer|normal)$")
			names(samples) = samples
			cancers = samples[grepl("cancer", samples)]
			normals = samples[grepl("normal", samples)]
			#cancers = "D02T_cancer"
		}
	} else {
		cancers = as.vector(can[,1])
		#print(cancers)
	}

	oneside.ram=F; min.out.conf=5;

	if (ctype %in% c("hg19")) ref="hg19"

	if (type == "germline") {

		x = lapply(samples, function(s) {
			######
			# skip the problematic samples
			if (s %in% c("kirc4856_normal", "aml2966_normal", "lv180_cancer", "lv180_normal")) return(NULL)

			out.fname = sprintf("%s/%s/cluster_%sm/%s.germline", dir, s, rasym, s) 
			if (file.exists(out.fname)) {
				print(paste("germline calls done already:", out.fname))
				return(NULL)
			}

			cldir = sprintf("%s/%s/cluster_%sm", dir, s, rasym)
			clfiles = list.files(cldir, "*.cluster$") 

			if ((rasym == "ra" || rasym == "um") && length(clfiles) != 24) {
				print(paste("skipping", s))
				return(NULL)
			} else if (rasym == "va" && length(clfiles) == 0) {
				print(paste("skipping", s))
				return(NULL)
			}

			if (rasym == "um" || rasym == "va") {
				oneside.ram = T
				min.out.conf = 2
			} else if (rasym == "ra") {
				oneside.ram = F
				min.out.conf = 5
			}

			print(paste("calling germline events for", s))
			y = call.germline(dir, s, ref=ref, rasym=rasym, oneside.ram=oneside.ram, min.out.conf=min.out.conf, mark.exo=F, contig=F, verbose=F)

			return(1)
		})

	} else if (type == "somatic") {

		x = lapply(samples, function(s) {

			if (s %in% cancers) {
				######
				# skip the problematic samples
				if (s %in% c("kirc4856_cancer", "aml2966_cancer")) return(NULL)

				out.fname = sprintf("%s/%s/cluster_%sm/%s.somatic.matched", dir, s, rasym, s) 
				if (file.exists(out.fname)) {
					print(paste("somatic call file already exists:", out.fname))
					return(NULL)
				}

				cldir = sprintf("%s/%s/cluster_%sm", dir, s, rasym)
				clfiles = list.files(cldir, "*.cluster$") 
		
				if (ctype == "brca_trip_negative") {
					matched.control = "normal" 
				} else {
					matched.control = as.vector(mat[,1])
					#print(matched.control)
				}

				cldir = sprintf("%s/%s/cluster_%sm", dir, matched.control, rasym)
				clfiles = list.files(cldir, "*.cluster.raw$") 

				print(paste("calling somatic events for", s))

				if (rasym != "ra") {
					nonmatched.controls = NULL
				} else {
					if (ref == "hg18") { 
						nonmatched.controls = nonmatched.controls.hg18 
					} else if (ref == "hg19") {
						nonmatched.controls = nonmatched.controls.hg19
						nonmatched.controls = nonmatched.controls[grep(matched.control, nonmatched.controls, invert=T)]
						#write.msg(paste("nonmatched files:: \"", nonmatched.controls, "\""))
					}
					names(nonmatched.controls) = nonmatched.controls
				}

				if (rasym == "um" || rasym == "va") {
					oneside.ram = T
					min.out.conf = 2
				} else if (rasym == "ra" || rasym == "youngte") {
					oneside.ram = F
					min.out.conf = 5
				}

				y = call.somatic(dir, s, rasym, matched.control, nonmatched.controls, ref=ref, oneside.ram=oneside.ram, min.out.conf=min.out.conf, verbose=F)
			}
		})
	}
}
