###############################################################
# Functions
###############################################################
# get.calls <- function(ctype, rasym, type, samples=NULL)
# get.somatic.batch <- function(dir, samples, rasym)
# get.somatic.va.batch <- function()
# get.va.somatic <- function(dir, s, rasym)
#
# analyze.va.somatic.exo
#
# region.overlap <- function(a, b, margin=0, single.chr=F)
# region.overlap.chr <- function(a, b, margin)
#
# make.igvplot.script <- function(df, session.fname, script.fname)
#  : to do
# comp.bgi.va <- function()
# comp.bgi.um <- function()
# load.bgi.calls <- function(dir, samples)
# make.igvsession.file <- function(ref="hg18", sample)
# 	: todo
###############################################################

source("/bitools/TEA/Tea-0.6.2/R/rid.r")

can <- read.table("cancer", header=F)
mat <- read.table("matched", header=F)

nonmatched.controls.hg18 = c("gbm0145_normal", "gbm0152_normal", "gbm0155_normal", "gbm0185_normal", "gbm0188_normal", "gbm0208_normal", "gbm0214_normal", "gbm0648_normal", "gbm0786_normal", "gbm0877_normal", "gbm0881_normal", "gbm1086_normal", "gbm1401_normal", "gbm1438_normal", "gbm1454_normal", "gbm1459_normal", "ov0723_normal", "ov0725_normal", "ov0751_normal", "ov0890_normal", "ov0980_normal", "ov1319_normal", "ov1411_normal", "na18506", "na18508")

if(file.exists("nonmatched")){
	nor <- read.table("nonmatched", header=F)
	nonmatched.controls.hg19 = as.vector(nor[,1])
} else{
	nonmatched.controls.hg19 = c()
}


library(RColorBrewer)

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
		print(cancers)
	}

	oneside.ram=F; min.out.conf=5; ref="hg18"

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
				print(paste(cldir))
				clfiles = list.files(cldir, "*.cluster$") 
				print(paste(clfiles))
		
				if (ctype == "brca_trip_negative") {
					matched.control = "normal" 
				} else {
					matched.control = as.vector(mat[,1])
					print(matched.control)
				}

				cldir = sprintf("%s/%s/cluster_%sm", dir, matched.control, rasym)
				clfiles = list.files(cldir, "*.cluster.raw$") 
				print(paste(clfiles))

				print(paste("calling somatic events for", s))

				if (rasym != "ra") {
					nonmatched.controls = NULL
				} else {
					if (ref == "hg18") { 
						nonmatched.controls = nonmatched.controls.hg18 
					} else if (ref == "hg19") {
						nonmatched.controls = nonmatched.controls.hg19
						nonmatched.controls = nonmatched.controls[grep(matched.control, nonmatched.controls, invert=T)]
						write.msg(paste("nonmatched files:: \"", nonmatched.controls, "\""))
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

brca <- function() {
	ph = read.delim("brca/phenotype/TCGA_brca_samples.txt", sep="\t", skip=1, header=T)
}

analyze.va.somatic.exo <- function()
{
	dir = "/files/CBMI/parklab/alee/ra2" # gbm/ov/pr/mm
  	samples = list.files(dir, "*_(cancer|normal)$")
    names(samples) = samples

	x = do.call(rbind, lapply(samples, function(s) {
		fname = sprintf("%s/%s/cluster_vam/%s.somatic.exo", dir, s, s) 
		if (file.exists(fname)) {
			sprintf("loading %s... ", fname)
			cl = read.delim(fname, sep="\t", as.is=T)
			sprintf("done\n")
			return(cl)
		} else {
			sprintf("no somatic.exo file %s exists", fname)
			return(NULL)
		}
	}))
	rownames(x) = NULL
	x = x[order(-x$ram),]

	bdir = "/groups/park/alee/vrs/tea_run"
	#ctypes = c("aml", "brca_trip_negative", "kirc", "liver", "luad", "lusc", "ucec")
	ctypes = c("aml", "kirc", "liver", "luad", "lusc", "ucec")

	y = do.call(rbind, lapply(ctypes, function(ct) {
		print(paste("processing", ct))
		dir = sprintf("%s/%s", bdir, ct)
		samples = list.files(dir, "*_(cancer|normal)$")
	
		samples = samples[grep("kirc4856", samples, invert=T)]
		x = get.somatic.batch(dir, samples, "va")	
		return(do.call(rbind, x))
	}))

#		x = do.call(rbind, lapply(samples, function(s) {
#        	fname = sprintf("%s/%s/cluster_vam/%s.somatic.exo", dir, s, s)
#			if (file.exists(fname)) {
#           		sprintf("loading %s... ", fname)
#           		cl = read.delim(fname, sep="\t", as.is=T)
#           	 	sprintf("done\n")
#           	 	return(cl)
#        	} else {
#           	 sprintf("no somatic.exo file %s exists", fname)
#           	 return(NULL)
#        	}
#		}))
#		return(x)
#    }))
}

get.outsamples <- function(dir, ctype)
{
	outsamples = NULL

    if (grepl("aml", dir)) {
    	# for aml use 10 normal genome
        if (rasym == "um") outsamples = c("aml2966_cancer", "aml2966_normal")

        if (rasym == "ra") outsamples = c("aml2967_cancer", "aml2966_normal", "aml2997_cancer", "aml2985_cancer",
            "aml2986_cancer", "aml2989_cancer")
    } else 
	if (grepl("luad", dir)) {
		if (rasym == "um") outsamples = c("luad4422_cancer", "luad4422_normal", "luad4510_cancer", "luad4510_normal")
	} else
	if (grepl("ucec", dir)) outsamples = c("ucecA054_cancer") else
    if (grepl("brca", dir)) outsamples = c("brcaA04Q_cancer", "brcaA0IJ_cancer") else
    if (grepl("lusc", dir)) {
        if (rasym == "ra") outsamples = c("lusc2596_cancer", "lusc2756_cancer")
    } else
    if (!is.null(ctype) && ctype == "cr") outsamples = c("cra00r_normal") else
    if (!is.null(ctype) && ctype == "pr") outsamples = c("pr0508_cancer", "pr0508_normal", "pr1701_cancer", "pr1701_normal", "pr1783_cancer") else
    if (!is.null(ctype) && ctype == "gbm") {
        # use all gbm and ov normal for nonmatched filtering not limited to 10 normals due to too many candidates.
        outsamples = c("gbm0145_cancer", "gbm0145_normal", "gbm0155_cancer", "gbm0155_normal", "gbm0188_cancer", "gbm0188_normal",
            "gbm0881_cancer", "gbm1086_cancer", "gbm1086_normal", "gbm1401_cancer", "gbm1454_cancer")
    } else
    if (!is.null(ctype) && ctype == "ov") {
        outsamples = c("ov0725_cancer", "ov0751_cancer", "ov0980_cancer", "ov0980_normal", "ov0982_cancer", "ov0982_normal")
    } 

	return(outsamples)
}

get.somatic.batch <- function(dir, rasym, ref, nm.controls=10, ctype=NULL)
{
	samples = dir(dir, "cancer$")
	nonmatched.controls = dir(dir, "normal$")

	# for gbm, ov, pr, mm, cr 
	if (!is.null(ctype)) {
		samples = samples[grep(paste("^", ctype, sep=""), samples)]
		#if (ctype == "gbm" || ctype == "ov") {
		#	nonmatched.controls = nonmatched.controls[grep("gbm|ov", nonmatched.controls)]
		#} else {
			nonmatched.controls = nonmatched.controls[grep(paste("^", ctype, sep=""), nonmatched.controls)]
		#}
	}

	outsamples = get.outsamples (dir, ctype) 

	if (!is.null(outsamples)) { 
		samples = setdiff(samples, outsamples)
		nonmatched.controls = setdiff(nonmatched.controls, outsamples)
	}

   	cll = lapply(samples, function(s) {
		print(paste("processing", s))
		out.fname = sprintf("%s/%s/cluster_%sm/%s.somatic", dir, s, rasym, s)

		if (!file.exists(out.fname)) {
   			if (grepl("cancer", s)) matched.control = sub("cancer", "normal", s) else 
			if (grepl("normal", s)) matched.control = sub("normal", "cancer", s) else
				stop(pasate("no matched control for sample", s))

			nonmatched.controls = setdiff(nonmatched.controls, matched.control)
			#if (length(nonmatched.controls) > nm.controls) nonmatched.controls = nonmatched.controls[1:nm.controls]

			print(paste("nonmatched controls:", paste(nonmatched.controls, collapse=",")))

			if (rasym == "va") cl = get.va.somatic(dir, s, mark.exo) else 
			if (rasym == "um") {
					oneside.ram = T
					min.acr = 0
					min.ram = 6
					min.acrr = 0.5
					min.out.conf = 3
					seqid=F;contig=T
					#call.somatic.lsf(dir, s, rasym, matched.control, nonmatched.controls, ref, oneside.ram=T, min.ram=min.ram,
                    #    min.acr=min.acr, min.acrr=min.acrr, matched.cram=1, matched.cacr=1, contig=T, seqid=T, min.out.conf=min.out.conf)
					cl = call.somatic(dir, s, rasym, matched.control, nonmatched.controls, ref, oneside.ram=oneside.ram, min.ram=min.ram, 
						min.acr=min.acr, min.acrr=min.acrr, matched.cram=1, matched.cacr=1, contig=contig, seqid=seqid, min.out.conf=min.out.conf)
			} else 
			if (rasym == "ra") {
					cl = call.somatic(dir, s, rasym, matched.control, nonmatched.controls, ref)
			}
		} else {
			cl = read.delim(out.fname, sep="\t", as.is=T)
		}
		return(cl)
	})
	return(cll)
}


get.va.somatic <- function(dir, s, mark.exo=T)
{
	chrl = load.ref.info("hg18", out.chrl=T)$chrl

	if (grepl("cancer", s)) matched.control = sub("cancer", "normal", s) else matched.control = sub("normal", "cancer", s)
    in.fname = sprintf("%s/%s/cluster_%sm/%s.cluster", dir, s, rasym, s)
    ram.fname = sprintf("%s/%s/cluster_ram/%s.cl.raw.RData", dir, s, s)

    control.fname = sprintf("%s/%s/cluster_%sm/%s.cluster.raw", dir, matched.control, rasym, matched.control)
    control.ram.fname = sprintf("%s/%s/cluster_ram/%s.cl.raw.RData", dir, matched.control, matched.control)

    out.fname = sprintf("%s/%s/cluster_%sm/%s.somatic", dir, s, rasym, s)
    out2.fname = sprintf("%s/%s/cluster_%sm/%s.somatic.exo", dir, s, rasym, s)

	if (file.exists(in.fname)) { 
		cl = read.delim(in.fname, sep="\t", header=T) 
	} else {
		sprintf("no cluster file: %s\n", in.fname)
		return(NULL)
	}

	if (nrow(cl) == 0) return(NULL)

	if (file.exists(control.fname)) {
		raw = read.delim(control.fname, sep="\t", header=T)
	} else {
		sprintf("no cluster file: %s\n", control.fname)
		return(NULL)
	}

	cl = cbind(sample=s, cl)
	control = region.overlap(cl, raw)
	cl = cl[control==0,]
	
	if (mark.exo && nrow(cl)>0) {
		x  = count.ram(cl, dir, s, chrl, fr.ram=0.2, max.ram=5, margin=0)
		if (!is.null(x)) {
			cl = x
			cl$exo = ifelse(cl$endo.ram <= cl$ram.cutoff, 1, 0)
		}
	}
	write.table(cl, out.fname, sep="\t", quote=F, row.names=F)
	if ("exo" %in% colnames(cl)) write.table(cl[cl$exo==1,], out2.fname, sep="\t", quote=F, row.names=F)
	
	return(cl)
}

region.overlap <- function(a, b, margin=0, single.chr=F)
{
    if (single.chr)  return(region.overlap.chr(a, b, margin))
    a$idx = 1:nrow(a)

    chrl = unique(as.character(a$chr))
    ov = do.call(rbind, lapply(chrl, function(ch) {
      print(paste("processing", ch))
      x = data.frame(idx=a$idx[a$chr==ch], ov=region.overlap.chr(a[a$chr==ch,], b[b$chr==ch,], margin))
      return(x)
    }))
    ov = ov[order(ov$idx),]
    return(ov$ov)
}

region.overlap.chr <- function(a, b, margin)
{
      ocnt = countOverlaps(IRanges(a$s-margin, a$e+margin), IRanges(b$s, b$e))
      ov = ifelse(ocnt>0, 1, 0)
      return(ov)
}

make.igvplot.script <- function(df, session.fname, script.fname)
{
	# make a script for plotting the regions in igv 

	# load session if available, otherwise, generate one
	if (!file.exists(session.fname)) {
		make.igvsession.file(session.fname)
	}
	cmd = sprintf("load %s", session.fname) 
}


comp.bgi.va <- function()
{
	dir = "/groups/park/alee/vrs/tea_run/liver"
  	samples = list.files(dir, "*_(cancer|normal)$")
    names(samples) = samples

	bgi = load.bgi.calls(dir, samples)

	teal = get.somatic.batch(dir, samples, "va") 
	tea = do.call(rbind, teal)
	tea.exo = do.call(rbind, lapply(teal, function(x) x[x$exo==1,]))
	rownames(tea) = rownames(tea.exo) = NULL

	# category: one/two side ram + one/two side cr
	tea$cat[tea$pram == 0 | tea$nram == 0] = 1 # oneside ram
	tea$cat[tea$pram > 0 & tea$nram > 0] = 2 # oneside ram

	tea.hbv = tea[tea$rep.repeat == 21326584,]
	tea.exo.hbv = tea.exo[tea.exo$rep.repeat == 21326584,]

	ov = region.overlap(tea.hbv, bgi)
	# number of calls depending to the ram and cr cutoff
	x = t(do.call(rbind, lapply(c(seq(2,5), seq(6,20,2)), function(n) {
		return(c(length(which(tea.hbv$ram>=n)), length(which(tea.hbv$ram>=n & tea.hbv$acr>0)), length(which(tea.hbv$ram>=n & tea.hbv$acr>1))))
	})))
	colnames(x) = c(seq(2,5), seq(6,20,2))	
	pdf("figure/HBV.events.freq.pdf")
	par(mgp=c(3, 1, 0))
	barplot(x, beside=T, xlab="Number of supporting paired reads", ylab="Number of events", col=three.col,border=F)
 	legend("topright", legend=c("#clipped reads >= 0", "#clipped reads >=1", "#clipped reads >=2"), fill=three.col, bty="n")
	dev.off()

	# category plot
	pdf("figure/category.ram2.pdf")
  	b = barplot(table(tea.hbv$cat), col="grey", space=0.5, names.arg=c("One side RAM", "Two side RAM"), ylab="Number of events") 
	dev.off()

	pdf("figure/category.ram3.cr.pdf")
  	b = barplot(table(tea.hbv$cat[tea.hbv$ram >=3 & tea.hbv$acr >0]), col="grey", space=0.5, names.arg=c("One side RAM", "Two side RAM"), ylab="Number of events") 
	dev.off()


	margin = 500
	x = do.call(rbind, lapply(samples, function(s) { 
		a = bgi[bgi$Sample.Library == s, ]
		b =	tea.hbv[tea.hbv$sample == s, ] 

		y = do.call(rbind, lapply(unique(as.character(a$chr)), function(chr) {
			a = a[a$chr == chr, ]
			b =	b[b$chr == chr,] 
			m = findOverlaps(IRanges(a$s, a$e) + margin, IRanges(b$s, b$e))@matchMatrix

			cnt = do.call(rbind, lapply(unique(m[,1]), function(i) {
				x = a$No..of.PE.read.support[i]
				idx = m[m[,1] == i,2]
				y = max(b$ram[idx])
				maxidx = which(b$ram[idx] == y)[1]
				z = b$ram[idx[maxidx]]+ b$acr[idx[maxidx]]
				print(paste(s, a$chr[i], a$s[i], a$e[i], b$chr[idx[maxidx]], b$s[idx[maxidx]], b$e[idx[maxidx]], x,y,z))
				return(c(x,y,z))
			}))
		return(cnt)
		}))
		return(y)
	}))
	pdf("figure/HBV.bgi.tea.ram.pdf")
	plot(x[,1], x[,2], xlab="BGI", ylab="Tea", main="Numer of supporting paired reads")
	abline(0, 1, lty=3)
	dev.off()
	
	# venn diagram with ram2 cutoff
	# need to do this againg considering the sample id; currently doing a pooled analysis: 
	# but all the insetions are unique across individuals so the current analysis is fine for now
	x = bgi[, c("chr", "s", "e")]
	y = tea.hbv[, c("chr", "s", "e")]
	m = rbind(x, y)
	um = get.usite(m)

	xv = which(region.overlap(um, x, margin=500)==1)
	yv = which(region.overlap(um, y, margin=500)==1)
 
	d = list(xv, yv)
	weights = draw.venn(d,  c("BGI", "Tea"), "figure/HBV.bgi.tea.margin500.ram2.pdf")
	
	# venn diagram with ram3 cutoff
	x = bgi[bgi$No..of.PE.read.support >=3,c("chr", "s", "e")]
	y = tea.hbv[tea.hbv$ram >=3, c("chr", "s", "e")]
	m = rbind(x, y)
	um = get.usite(m)

	xv = which(region.overlap(um, x, margin=500)==1)
	yv = which(region.overlap(um, y, margin=500)==1)

	d = list(xv, yv)
	weights = draw.venn(d,  c("BGI", "Tea"), "figure/HBV.bgi.tea.margin500.ram3.pdf")

  	x = bgi[, c("chr", "s", "e")]
    y = tea.hbv[tea.hbv$acr>0, c("chr", "s", "e")]
    m = rbind(x, y)
    um = get.usite(m)

    xv = which(region.overlap(um, x, margin=500)==1)
    yv = which(region.overlap(um, y, margin=500)==1)
	d = list(xv, yv)

	weights = draw.venn(d,  c("BGI", "Tea"), "figure/HBV.bgi.tea.margin500.ram2.cr.pdf")

  	x = bgi[bgi$No..of.PE.read.support >=3, c("chr", "s", "e")]
    y = tea.hbv[tea.hbv$acr>0, c("chr", "s", "e")]
    m = rbind(x, y)
    um = get.usite(m)

    xv = which(region.overlap(um, x, margin=500)==1)
    yv = which(region.overlap(um, y, margin=500)==1)
	d = list(xv, yv)

	weights = draw.venn(d,  c("BGI", "Tea"), "figure/HBV.bgi.tea.margin500.ram3.cr.pdf")

}

comp.bgi.um <- function()
{
    dir = "/groups/park/alee/vrs/tea_run/liver"
    samples = list.files(dir, "*_(cancer|normal)$")
    names(samples) = samples

	samples = samples[-grep("lv70", samples)]

    bgi = load.bgi.calls(dir, samples)
    teal = get.somatic.batch(dir, samples, rasym="um")
    tea = do.call(rbind, teal)
    rownames(tea) =  NULL

    tea$cat[tea$pram == 0 | tea$nram == 0] = 1 # oneside ram
    tea$cat[tea$pram > 0 & tea$nram > 0] = 2 # oneside ram

	ov = region.overlap(bgi, tea, margin=500)
	tea$ov = region.overlap(tea, bgi, margin=500)
	tea.ordered = tea[order(-tea$ram),]
    # category: one/two side ram + one/two side cr
    margin = 500
    cnt = do.call(rbind, lapply(samples, function(s) {
        a = bgi[bgi$Sample.Library == s, ]
        b = tea[tea$sample == s, ]

        y = do.call(rbind, lapply(unique(as.character(a$chr)), function(chr) {
            a = a[a$chr == chr, ]
            b = b[b$chr == chr,]
            m = findOverlaps(IRanges(a$s, a$e) + margin, IRanges(b$s, b$e))@matchMatrix

            cnt = do.call(rbind, lapply(unique(m[,1]), function(i) {
                x = a$No..of.PE.read.support[i]
                idx = m[m[,1] == i,2]
                y = max(b$ram[idx])
                maxidx = which(b$ram[idx] == y)[1]
                z = b$acr[idx[maxidx]]
                print(paste(s, a$chr[i], a$s[i], a$e[i], b$chr[idx[maxidx]], b$s[idx[maxidx]], b$e[idx[maxidx]], x,y,z))
                return(c(x,y,z))
            }))
        return(cnt)
        }))
        return(y)
    }))
    pdf("figure/um.bgi.tea.ram.pdf")
    plot(cnt[,1], cnt[,2], xlab="BGI", ylab="Tea", main="Numer of supporting paired reads")
    abline(0, 1, lty=3)
    dev.off()

    # venn diagram with ram2 cutoff
    # need to do this againg considering the sample id; currently doing a pooled analysis: 
    # but all the insetions are unique across individuals so the current analysis is fine for now
    x = bgi[bgi$No..of.PE.read.support>=6, c("chr", "s", "e")]
    y = tea[tea$ram >= 10 & tea$acr >=3, c("chr", "s", "e")]
    m = rbind(x, y)
    um = get.usite(m)

    xv = which(region.overlap(um, x, margin=500)==1)
    yv = which(region.overlap(um, y, margin=500)==1)

    d = list(xv, yv)
    weights = draw.venn(d,  c("BGI", "Tea"), "figure/um.bgi.tea.margin500.ram6.acr3.pdf")

}

get.usite <- function(df)
{
  us = data.frame(do.call(rbind, lapply(unique(df$chr), function(chr) {
    if (nrow(df[df$chr==chr,]) >0) {
      r = reduce(IRanges(df$s[df$chr == chr], df$e[df$chr == chr]))
      return(data.frame(chr=chr, s=r@start, e=r@start+r@width-1))
    } else return(NULL)
  })))
    rownames(us) = NULL
  return(us)
}

draw.venn <- function(d, setNames, fname) {
    require(Vennerable)
    pdf(fname)
    venn = Venn(d, SetNames=setNames)
    weights = sapply(venn@IntersectionSets, length)
    Vcombo<-compute.Venn(venn)
   
    vennSchemaColori<-VennThemes(Vcombo, colourAlgorithm="sequential")
   
    # intersectino regions colors
    vennSchemaColori$Face$'001'$fill<-"#fe8081"
    vennSchemaColori$Face$'010'$fill<-"#c2d2a3"
    vennSchemaColori$Face$'100'$fill<-"#99bae5"
    vennSchemaColori$Face$'101'$fill<-"#987ba7"
    vennSchemaColori$Face$'011'$fill<-"#c29264"
    vennSchemaColori$Face$'110'$fill<-"#7aa3b5"
    vennSchemaColori$Face$'111'$fill<-"#7e849a"
   
    # circles border color
    vennSchemaColori$Set$Set1$col<-"black"
    vennSchemaColori$Set$Set2$col<-"black"
    vennSchemaColori$Set$Set3$col<-"black"
   
    # label colors
    vennSchemaColori$SetText$Set1$col<-"black"
    vennSchemaColori$SetText$Set2$col<-"black"
    vennSchemaColori$SetText$Set3$col<-"black"

    plot(Vcombo, gpList=vennSchemaColori)
    dev.off()
    return(weights)
}


load.bgi.calls <- function(dir, samples)
{

	bgi = read.delim(sprintf("%s/ng.2295-S3.txt", dir), sep="\t", as.is=T)
	bgi = bgi[bgi$hg18.by.Liftover != "",]
	bgi$Sample.Library = sub("^", "lv", bgi$Sample.Library)
	bgi$Sample.Library = sub("T", "_cancer", bgi$Sample.Library)
	bgi$Sample.Library = sub("N", "_normal", bgi$Sample.Library)
	bgi = bgi[bgi$Sample.Library %in% samples,]
	x = data.frame(do.call(rbind, strsplit(bgi$hg18.by.Liftover, "[:-]")), stringsAsFactors=F)
	colnames(x) = c("chr", "s", "e")
	bgi = cbind(x, bgi)
	bgi$s = as.integer(bgi$s)
	bgi$e = as.integer(bgi$e)
	return(bgi)
}

make.igvsession.file <- function(ref="hg18", sample)
{
#	dir = "/home/el114/repeat_analysis/igv/session"
#	templates = list(hg18="a", hg19="b")(	)

#	fname = sprintf("%s/...", dir) 
#	cmd = "cp 
#	system("perl "

}
