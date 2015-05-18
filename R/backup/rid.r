#
# Alice E. Lee: ejalice.lee@gmail.com
#
#################################################################################################
# changes
#################################################################################################
#
# rid
#		- reference support: hg18/hg19
#		- endogeneoues/exogenous
#		- allow one side ram 
#		- assembly prefix (ra, va, etc.)
#      	- cbam.chr to handle the referecnce without "chr"
#       - rmasker.filter.margin
#       - merge.family: merge TEs of the same family in cluster detection
#       - bug fixed:
#           clipping read position considering base insertions (cigar string I)
#			a bug to omit the clusters when the pairing does not meet the cluster size criteria
#
# rid.te.primate 
#		- reference support
#			* hg18/hg19 (human), ponAbe2 (orangutan), panTro3 (chimpanzee), rheMac2 (Rhesus)  
#		- merge.family 
#			* te library for each te family along with family-wise filtering of known instances
# 
# rid.te (original; need to check the integrity due to some further changes)
#
# call.somatic
#		- add a filtering based on the fraction of reads instead of the read counts
#
# s2n 
#       - calculate signal to noise ratio 
#       - 2*sqrt(#pram_downstream_of_nbp * #nram_upstream_of_pbp) - (#pram_upstreadm_of_nbp + #nram_downstream_of_pbp)
#################################################################################################

library(IRanges)
library(Rsamtools)
library(spp)

tea.base = Sys.getenv("tea_base")
if (tea.base == "") stop("tea_base environment variable needs to be set")
sprintf("tea.base: %s", tea.base)

# rid can detect both endogeneous and exogeneous insertions 
# rasym needs to be set up to represent different sequence libraries 
# (example)
# ra: TE sequence libraries
# va: virus sequence libraries
# um: no sequence library was used but identifying the clusters of unmapped reads
rid <- function(sample, dir, chr=NULL, no.clipped=F, ref="hg18", rasym="ra", min.ram=3, oneside.ram=F, exo=F, cbam.chr=T, jittering=2, merged.rmasker=T, merge.family=T, rmasker.filter.margin=500, annot.oi=T, annot.gene=T, verbose=F)
{
	# input parameter integrity check
	ref.support = c("hg18", "hg19")
	if (! ref %in% ref.support)	stop(sprintf("The reference %s is not supported\n", ref))

	# need to orgnize the relation between merged.rmasker and merge.family	
	merged.rmasker = merge.family

	if (exo) {
		oneside.ram = T
		merge.family = F
		annot.oi = F
	}
	
	if (rasym == "um") stringent.pair=T else stringent.pair=F

	# load the annotations
	ref.annot = load.ref.info(ref, out.chrl=T, out.ril=!exo&annot.oi, out.gene=annot.gene, out.gap=F, merged.rmasker=merged.rmasker, verbose) 

	# input files
	bam.dir = sprintf("%s/%s/bam", dir, sample)

	ram.file = sprintf("%s/%s.%sm.bz2", bam.dir, sample, rasym)
	rl.file = sprintf("%s/%s.rl", bam.dir, sample)
	isize.file = sprintf("%s/%s.isize", bam.dir, sample)

   	if (no.clipped) {
       cbam.file = NULL
    } else {
		cbam.file = sprintf("%s/%s.sorted.softclips.consd.bam", bam.dir, sample)
		if (!file.exists(cbam.file)) cbam.file =  sprintf("%s/%s.softclips.consd.bam", bam.dir, sample)
		if (!file.exists(cbam.file)) stop("there is no clipped bam file")
	}

	# output files
	cl.dir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)
	if (!file.exists(cl.dir)) { dir.create(cl.dir, recursive=T, mode="0755")}

	if (is.null(chr)) {
		cl.prefix = sprintf("%s/%s", cl.dir, sample)
	} else {
	  cl.prefix = sprintf("%s/%s.%s", cl.dir, sample, chr)
	}
	cl.raw.file = sprintf("%s.cluster.raw", cl.prefix) 
	cl.file = sprintf("%s.cluster", cl.prefix) 
	clipped.file = sprintf("%s.clipped", cl.prefix)

	write.msg(paste("loading rl and isize from", rl.file, isize.file, "..."))
	rl = load.rl(rl.file)
	is = load.isize(isize.file, rl)
	if (is.null(rmasker.filter.margin)) rmasker.filter.margin = is$ins.margin
	write.msg(sprintf("fragment: %d (mu: %d, sd: %d), intra.gap: %d, inter.gap: %d, ins.margin: %d, rmasker.filter.margin: %d", is$fr, is$mu, is$sd, is$intra.gap, is$inter.gap, is$ins.margin, rmasker.filter.margin))

	write.msg(paste("loading rams from", ram.file, "..."))
	ram = load.ram(ram.file, separated=F)
	write.msg("done loading rams.")

	if (!is.null(chr)  && is.null(ram[[chr]]))  { stop(sprintf("there is no ram in %s", chr)) }

	if (is.null(chr))   chrl = ref.annot$chrl else chrl = c(chr)
	names(chrl) = chrl

	# write the headers
	colh = get.col.header()
	write(paste(colh$cl.raw, collapse="\t"), cl.raw.file)
	write(paste(colh$cl, collapse="\t"), cl.file)
	write(paste(colh$clipped, collapse="\t"), clipped.file)
	
	ptm0 = proc.time()
	cl = lapply(chrl, function(c) {
		write.msg(paste("processing", c))
		ptm = proc.time()
	  	if (is.null(ram[[c]]))	return(NULL)
		if (!oneside.ram && ( dim(ram[[c]]$p)[1] == 0 || dim(ram[[c]]$m)[1] == 0)) return (NULL)

		p.cl = NULL 
		if ( nrow(ram[[c]]$p) > 0) {
		  	write.msg(paste("extracting positive strand ram (pram) clusters ..."))
		  	ram[[c]]$p$repeatname = gsub("/", "_", ram[[c]]$p$repeatname) # replace ALR/ALPHA to ALR_ALPHA
		 	p.cl = get.cluster2(ram[[c]]$p, 1, is$intra.gap, merge.family, verbose)
		  	write.msg(paste("done: positive strand clusters:", dim(p.cl)[1]))
			org.header=colnames(p.cl)
		}

		m.cl = NULL
		if ( nrow(ram[[c]]$m) > 0) {
		  	write.msg(paste("extracting negative strand  ram (nram) clusters ..."))
		  	ram[[c]]$m$repeatname = gsub("/", "_", ram[[c]]$m$repeatname) # replace ALR/ALPHA to ALR_ALPHA
		  	m.cl = get.cluster2(ram[[c]]$m, -1, is$intra.gap, merge.family, verbose)
		  	write.msg(paste("done: negative strand clusters:", dim(m.cl)[1]))
			org.header=colnames(m.cl)
		}

		pm.cl = NULL; paired.pidx = paired.midx = NULL 
		if (!is.null(p.cl) && !is.null(m.cl)) {
		  	write.msg(paste("pairing clusters ..."))
		  	pov = pair.cluster3(p.cl, m.cl, is$inter.gap, stringent.pair, verbose)
		  	pm.cl = pov$pcl
			paired.pidx = pov$ov[,1]
			paired.midx = pov$ov[,2]
			write.msg(paste("done: paired clusters:", dim(pm.cl)[1]))
		}

		if (!is.null(pm.cl) && nrow(pm.cl) > 0)  {
	    	pm.cl$s = pm.cl$s1
	    	pm.cl$e = pm.cl$e2+rl
		
	    	# define ram cluster boundaries and remove clusters whose boundary size is smaller than 2 * rl + 10 
	    	idx = which((pm.cl$e - pm.cl$s + 1) >= (2 * rl + 10))
	    	pm.cl = pm.cl[idx,]
			paired.pidx = pov$ov[idx,1]
			paired.midx = pov$ov[idx,2]
		}

		if (!is.null(p.cl)) colnames(p.cl) = paste(colnames(p.cl), '1', sep="")
	    if (!is.null(m.cl)) colnames(m.cl) = paste(colnames(m.cl), '2', sep="")
	
		ponly = NULL	
		if (!is.null(p.cl)) {	
		    empty.len = length(setdiff(1:dim(p.cl)[1], paired.pidx))
		    empty.num = rep(0, empty.len)
		    empty.chr = rep("", empty.len)
		    empty.df  = data.frame(empty.num, empty.num, empty.num, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, stringsAsFactors=F)
		    colnames(empty.df) = paste(org.header, '2', sep="")
	
		    ponly = cbind(p.cl[setdiff(1:dim(p.cl)[1], paired.pidx),], empty.df)
		    ponly$s = ponly$s1
		    ponly$e = ponly$e1 + is$fr
		    write.msg(paste("positive ram only clusters:", dim(ponly)[1]))
		}	

		monly = NULL	
		if (!is.null(m.cl)) {
		    empty.len = length(setdiff(1:dim(m.cl)[1], paired.midx))
		    empty.num = rep(0, empty.len)
		    empty.chr = rep("", empty.len)
		    empty.df  = data.frame(empty.num, empty.num, empty.num, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, stringsAsFactors=F)
		    colnames(empty.df) = paste(org.header, '1', sep="")
	
		    monly = cbind(empty.df, m.cl[setdiff(1:dim(m.cl)[1], paired.midx),])
		    monly$s = monly$s2 - is$fr
		    monly$e = monly$e2
		    write.msg(paste("negative ram only clusters:", dim(monly)[1]))
		}
	
	    cl.ch = rbind(pm.cl, ponly, monly)
	    cl.ch$chr = c
	    cl.ch$size = cl.ch$e - cl.ch$s + 1
	    cl.ch$ram = cl.ch$ram1 + cl.ch$ram2
	    cl.ch$rep.repeat = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	       return(paste(unique(c(strsplit(as.character(cl.ch$rep.repeat1[i]), ",")[[1]], strsplit(as.character(cl.ch$rep.repeat2[i]), ",")[[1]])), collapse=","))
	    }))
	    cl.ch$family = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	      if (cl.ch$family1[i] == "") return(cl.ch$family2[i]) else if
	         (cl.ch$family2[i] == "") return(cl.ch$family1[i]) else
	        return(cl.ch$family1[i]) # otherwise family1 is equal to family2  
	    }))
	    cl.ch$class = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	      if (cl.ch$class1[i] == "") return(cl.ch$class2[i]) else if
	         (cl.ch$class2[i] == "") return(cl.ch$class1[i]) else
	        return(cl.ch$class1[i]) # otherwise family1 is equal to family2 
	    }))
	
		cl.ch = cl.ch[, colh$cl.raw]

	    write.table(cl.ch, cl.raw.file, sep="\t", row.names=F, col.names=F, quote=F, append=T)
	    write.msg(sprintf("done writing %d raw clusters for %s", nrow(cl.ch),  c))
	
		colnames(cl.ch)[9:14] = c("pram", "nram", "pram_start", "pram_end", "nram_start", "nram_end")

		# remove clusters whose boundary size is smaller than 2 * rl + 10 ; related with update on the count.clipped2 (-rl)
	    cl.ch = cl.ch[cl.ch$size >= (2 * rl + 10), ]

		# filter out the clusters with less than minimum rams
   		cl.ch = cl.ch[cl.ch$ram >= min.ram, ]
		if (nrow(cl.ch) == 0 ) {
			write.msg(sprintf("no cluster with %d rams", min.ram))
			return(NULL)
		}
		if (!oneside.ram) cl.ch = cl.ch[cl.ch$pram >= 1 & cl.ch$nram >= 1, ]

	    if (no.clipped) {
	      clipped.cnt = data.frame(acr=0, pacr=0, nacr=0, cr=0, pcr=0, ncr=0, pbp=-1, nbp=-1, tsd=-9999, arr=0)
	    } else {
	      write.msg("counting clipped reads ...")
	      clipped = count.clipped2(cl.ch, cbam.file, verbose, jittering, cbam.chr, rl=rl)
	      clipped.cnt = data.frame(do.call(rbind, lapply(clipped, function(x) unlist(x[["ccnt"]]))))
	      write.msg(paste("done:", sum(clipped.cnt$cr), "clipped reads"))
	    }
		# for compatibility, update the clipped.cnt column header 
		colnames(clipped.cnt) = c("acr", "pacr", "nacr", "cr", "pcr", "ncr", "pbp", "nbp", "tsd", "acrr")	
		cl.ch = cbind(cl.ch, clipped.cnt)

	    # calculate scores
	    cl.ch$score = get.score(cl.ch, 20, 10)

		cl.ch$oi = "-" # not annotated
		if (!exo && annot.oi) {
			write.msg("annotating known repeat instances ...")
			cl.ch = annotate.oi(cl.ch, ref.annot$ril, rmasker.filter.margin, verbose)
			write.msg("done annotating.")
		}

		cl.ch$pgene = cl.ch$ngene = "-" 
		if (annot.gene) {
			write.msg("annotating genes ...")
			cl.ch = annotate.genes(cl.ch, ref.annot$genes, margin=2)
			write.msg("done annotating.")
		}

		if (exo && rasym != "um") cl.ch$desc = annotate.virus(cl.ch) else cl.ch$desc = "-"
	    cl.ch = cl.ch[, colh$cl]

	    write.table(cl.ch, cl.file, sep="\t", quote=F, row.names=F, append=T, col.names=F)
	    write.msg(sprintf("done writing %d clusters for %s", nrow(cl.ch), c))

	    if (no.clipped) {
	      clipped.detail=NULL
	    } else {
	      clipped.detail = data.frame(do.call(rbind, lapply(clipped, function(x)
	       if (!is.null(x$cdf)) return(x[["cdf"]]))))
	      if (dim(clipped.detail)[1] >0) {
	        colnames(clipped.detail) = c("chr", "s", "e", "cpos", "cigar", "tname", "ref.seq", "clipped.seq", "clipped.qual", "aligned")
	        clipped.detail = clipped.detail[, colh$clipped]
	        write.table(clipped.detail, clipped.file,  row.names=F, quote=F, sep="\t", append=T, col.names=F)
	        write.msg(sprintf("done writing %d clipped sequences for %s", nrow(clipped.detail), c))
	      }
	    }
		rm(cl.ch, clipped.detail)
		
	    write.msg(paste(c, "elapsed:", (proc.time()-ptm)[["elapsed"]]))

	    return(1)
	})
	rm (cl, ram)
	write.msg(paste("all elapsed:", (proc.time()-ptm0)[["elapsed"]]))
	
	return(1)
}

merge.cluster <- function(cl.ch)
{
	families = unique(cl.ch$family); names(families) = families
	x = do.call(rbind, lapply(families, function(f) {
		cl.ch.f = cl.ch[cl.ch$family == f, , drop=FALSE]
		return(merge.cluster.family(cl.ch.f)) 
	}))
	x = x[order(as.numeric(rownames(x))),]
	return(x)
}

merge.cluster.family <- function(cl.ch.f) 
{
	# process pbp first
	bps = unique(cl.ch.f$pbp) 	
	bps = tail(bps, -1)	# ignore -1
	cl.ch.f = data.frame(do.call(rbind, lapply(bps, function(bp) { 
		x = cl.ch.f[cl.ch.f$pbp == bp, ,drop=FALSE]
		if (nrow(x) == 1) return(x)
		return(get.merged.cluster(x))
	})))

	bps = unique(cl.ch.f$nbp) 	
	bps = tail(bps, -1)	# ignore -1
	cl.ch.f = data.frame(do.call(rbind, lapply(bps, function(bp) { 
		x = cl.ch.f[cl.ch.f$nbp == bp, ,drop=FALSE]
		if (nrow(x) == 1) return(x)
		return(get.merged.cluster(x))
	})))
	return(cl.ch.f)
}

get.merged.cluster <- function(x) 
{
   # choose a cluster with the maxium nacr counts 
        i = which(x$nacr == max(x$nacr))[1]
        xx = x[i,]
        xx$s = min(x$s)
        xx$s = max(x$e) 
        xx$pram_start = min(x$pram_start)
        xx$pram_end = max(x$pram_end)
        xx$nram_start = min(x$nram_start)
        xx$nram_end = max(x$nram_end)

        xx$size = xx$e - xx$s + 1
        xx$tsd  = xx$nbp - xx$pbp - 1

        xx$ram = sum(x$ram)
        xx$pram = sum(x$pram)
        xx$nram = sum(x$nram)


        xx$ncr = x$ncr[i]
        xx$cr = xx$pcr + xx$ncr
        xx$nacr = x$nacr[i]
        xx$acr = xx$pacr + xx$nacr

        xx$acrr = round(xx$acr / xx$cr, 2)

        return(xx)
}

get.col.header <- function(annot.gap=F)
{
	colh = list()
	colh$cl.raw = c("chr", "s", "e", "size", "rep.repeat", "family", "class", "ram", "ram1", "ram2", "s1", "e1", 
			"s2", "e2", "pos1", "pos2", "rep.repeat1", "rep.repeat2", "repeat.name1", "repeat.name2", "rname1", "rname2")
	colh$cl = c("chr", "s", "e", "size", "tsd", "pbp", "nbp", "rep.repeat", "family", "class", "ram", "pram", "nram", 
			"cr", "pcr", "ncr", "acr", "pacr", "nacr", "acrr", "pram_start", "pram_end", "nram_start", "nram_end", "pgene", "ngene", "score", "oi", "desc")
	if (annot.gap) colh$cl = c(colh$cl, "gap")
	colh$clipped = c("chr", "s", "e", "cpos", "aligned", "cigar", "tname", "ref.seq", "clipped.seq", "clipped.qual")

	return(colh)
}

annotate.virus <- function(cl) 
{
	vrs.annot = load.virusannot()	
	vrs = unlist(lapply(1:nrow(cl), function(i) {
		vs = strsplit(cl$rep.repeat[i], "_&_")[[1]]		
		return(paste(sapply(vs, function(v) return(vrs.annot$name[vrs.annot$id == v])), collapse=","))
		
		vrs.annot$id == vs
	}))
	return(vrs)
}

load.virusannot <- function()
{
	annot.file = sprintf("%s/lib/viruses.txt", tea.base)

	colClasses = c("NULL", "character", "NULL", "NULL", "character")
    annot = read.table(annot.file, sep="|", header=F, quote="", colClasses = colClasses, as.is=T)
	colnames(annot) = c("id", "name")

	annot$name = sub("^\\s+", "", annot$name)
	annot$name = sub(", (complete |)genome", "", annot$name) 
	
    return(annot)
}

annotate.oi <- function(cl, ril, margin, verbose)
{
	write.msg(sprintf("marking clusters near known repeat instances: margin %d", margin))

    chrl = unique(as.character(cl$chr)); names(chrl) = chrl
    df = data.frame(do.call(rbind, lapply(chrl, function(ch) {
      if (verbose) print(paste("processing", ch))
		cl.ch = cl[cl$chr == ch,]
		cl.ch$oi = unlist(lapply(1:nrow(cl.ch), function(i) {
			if (verbose & i%%100 ==0) print(paste("processing", i))
			repeats = unlist(strsplit(strsplit(cl.ch$rep.repeat[i], "_&_")[[1]], ",")[[1]]) # added _&_ for the merged virus names
			# the annotation point was defined by pbp or nbp 
			# if neither pbp nor nbp was defined, mid point between s and e was used
			if (cl.ch$pbp[i] != -1 && cl.ch$nbp[i] != -1) {
				s = min(cl.ch$pbp[i], cl.ch$nbp[i])
				e = max(cl.ch$pbp[i], cl.ch$nbp[i])
			} else if (cl.ch$pbp[i] != -1 && cl.ch$nbp[i] == -1) {
				s = e = cl.ch$pbp[i]
			} else if (cl.ch$pbp[i] == -1 && cl.ch$nbp[i] != -1) {
				s = e = cl.ch$nbp[i]
			} else { # neither pbp nor nbp was defined
				s = e = floor((cl.ch$s[i] + cl.ch$e[i])/2)
			}
			oi = min(unlist(lapply(repeats, function(r) {
				oi = ifelse(countOverlaps(IRanges(s, e) + margin, IRanges(ril[[r]][[ch]]$s, ril[[r]][[ch]]$e))==0, 1, 0)
			})))
			return(oi)
		}))
		return(cl.ch)
	})))
	rownames(df) = NULL
	return(df)
}

annotate.oi.chr <- function(cl.ch, ril.ch, margin, verbose)
{
    write.msg(sprintf("marking clusters near known repeat instances: margin %d", margin))

    oi = unlist(lapply(1:nrow(cl.ch), function(i) {
		if (verbose & i%%100 ==0) print(paste("processing", i))
		repeats = unlist(strsplit(strsplit(cl.ch$rep.repeat[i], "_&_")[[1]], ",")[[1]]) # added _&_ for the merged virus names
		# the annotation point was defined by pbp or nbp 
		 # if neither pbp nor nbp was defined, mid point between s and e was used
            if (cl.ch$pbp[i] != -1 && cl.ch$nbp[i] != -1) {
                s = min(cl.ch$pbp[i], cl.ch$nbp[i])
                e = max(cl.ch$pbp[i], cl.ch$nbp[i])
            } else if (cl.ch$pbp[i] != -1 && cl.ch$nbp[i] == -1) {
                s = e = cl.ch$pbp[i]
            } else if (cl.ch$pbp[i] == -1 && cl.ch$nbp[i] != -1) {
                s = e = cl.ch$nbp[i]
            } else { # neither pbp nor nbp was defined
                s = e = floor((cl.ch$s[i] + cl.ch$e[i])/2)
            }
            oi = min(unlist(lapply(repeats, function(r) {
                oi = ifelse(countOverlaps(IRanges(s, e) + margin, IRanges(ril.ch[[r]]$s, ril.ch[[r]]$e))==0, 1, 0)
            })))
            return(oi)
	}))
    return(oi)
}


rid.te.primate <- function(sample, dir, chr=NULL, no.clipped=F, ref, rasym=NULL, min.ram=3, oneside.ram=F, exo=F, cbam.chr=T, jittering=2, merged.rmasker=F, merge.family=F, rmasker.filter.margin=500, annot.oi=T, annot.gene=F, annot.gap=T, verbose=F)
{
	# load the annotations
	ref.annot = load.ref.info(ref, out.chrl=T, out.ril=T, out.gene=annot.gene, out.gap=annot.gap, merged.rmasker=merged.rmasker, verbose) 
	ril = ref.annot$ril

	if (is.null(chr))   chrl = ref.annot$chrl else  chrl = c(chr)
	names(chrl) = chrl

	# input files: ram, clipped reads, isize info) 
	bam.dir = sprintf("%s/%s/bam", dir, sample)
	ram.file = sprintf("%s/%s.ram.bz2", bam.dir, sample)
	cbam.file = sprintf("%s/%s.sorted.softclips.consd.bam", bam.dir, sample)
	if (no.clipped) cbam.file = NULL
	rl.file = sprintf("%s/%s.rl", bam.dir, sample)
	isize.file = sprintf("%s/%s.isize", bam.dir, sample)

	# load rams
	write.msg(paste("loading rams from", ram.file, "..."))
	ram = load.ram(ram.file, separated=F)
	write.msg("done loading rams.\n")

	# load isize and read length
	write.msg(paste("loading rl and isize from", rl.file, isize.file, "..."))
	rl = load.rl(rl.file)
	is = load.isize(isize.file, rl)
	write.msg(paste("fragment:", is$fr, "(mu:", is$mu, "sd:", is$sd, ") intra.gap:", is$intra.gap, " inter.gap:", is$inter.gap, " ins.margin:", is$ins.margin))

	# ouput files: raw.cluser, cluster, clipped files 
	if (is.null(rasym)) cl.dir = sprintf("%s/%s/cluster", dir, sample) else
		cl.dir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)
	if (!file.exists(cl.dir)) dir.create(cl.dir, recursive=T, mode="0755")
	if (is.null(chr)) out.prefix = sprintf("%s/%s", cl.dir, sample) else
	  out.prefix = sprintf("%s/%s.%s", cl.dir, sample, chr)

	cl.raw.file = sprintf("%s.cluster.raw", out.prefix) 
	cl.file = sprintf("%s.cluster", out.prefix) 
	clipped.file = sprintf("%s.clipped", out.prefix)

	# write headers
	colh = get.col.header(annot.gap = annot.gap)
	write(paste(colh$cl.raw, collapse="\t"), cl.raw.file)
	write(paste(colh$cl, collapse="\t"), cl.file)
	write(paste(colh$clipped, collapse="\t"), clipped.file)

	ptm0 = proc.time()
	cl = lapply(chrl, function(c) {
		write.msg(paste("processing", c))
		ptm = proc.time()
		if (is.null(ram[[c]]))  return(NULL)
		if (!oneside.ram && ( dim(ram[[c]]$p)[1] == 0 || dim(ram[[c]]$m)[1] == 0)) return (NULL)
	 
		p.cl = NULL
		if (nrow(ram[[c]]$p) > 0) {
			write.msg(paste("extracting positive strand clusters ..."))
			ram[[c]]$p$repeatname = gsub("/", "_", ram[[c]]$p$repeatname) # (e.g.) ALR/ALPHA to ALR_ALPHA to prevent errors
		 	p.cl = get.cluster2(ram[[c]]$p, 1, is$intra.gap, merge.family, verbose)
			write.msg(paste("done: positive strand clusters:", dim(p.cl)[1]))
			org.header = colnames(p.cl)
		}

		m.cl = NULL
		if ( nrow(ram[[c]]$m) > 0) {
			write.msg(paste("extracting negative strand clusters ..."))
			ram[[c]]$m$repeatname = gsub("/", "_", ram[[c]]$m$repeatname) 
		 	m.cl = get.cluster2(ram[[c]]$m, 1, is$intra.gap, merge.family, verbose)
			write.msg(paste("done: negative strand clusters:", dim(m.cl)[1]))
		}

		pm.cl = NULL; paired.pidx = paired.midx = NULL
		 if ( !is.null(p.cl) && !is.null(m.cl)) {
			write.msg(paste("pairing clusters ..."))
			#pov = pair.cluster2(p.cl, m.cl, is$inter.gap)
			pov = pair.cluster3(p.cl, m.cl, is$inter.gap)
			pm.cl = pov$pcl
			write.msg(paste("done: paired clusters:", dim(pm.cl)[1]))
		}
	
		if (nrow(pm.cl) > 0)	{ 
		  	pm.cl$s = pm.cl$s1
		  	pm.cl$e = pm.cl$e2+rl
		}

		colnames(p.cl) = paste(colnames(p.cl), '1', sep="")
		colnames(m.cl) = paste(colnames(m.cl), '2', sep="")

		empty.len = length(setdiff(1:dim(p.cl)[1], pov$ov[,1]))
		empty.num = rep(0, empty.len)
		empty.chr = rep("", empty.len)
		empty.df  = data.frame(empty.num, empty.num, empty.num, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr)
		colnames(empty.df) = colnames(m.cl)

		ponly = cbind(p.cl[setdiff(1:dim(p.cl)[1], pov$ov[,1]),], empty.df)
		ponly$s = ponly$s1
		ponly$e = ponly$e1 + is$fr
		write.msg(paste("positive ram only clusters:", dim(ponly)[1]))

		empty.len = length(setdiff(1:dim(m.cl)[1], pov$ov[,2]))
		empty.num = rep(0,empty.len)
		empty.chr = rep("", empty.len)
		empty.df  = data.frame(empty.num, empty.num, empty.num, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr)
		colnames(empty.df) = colnames(p.cl)
		monly = cbind(empty.df, m.cl[setdiff(1:dim(m.cl)[1], pov$ov[,2]),])
		monly$s = monly$s2 - is$fr
		monly$e = monly$e2
		write.msg(paste("negative ram only clusters:", dim(monly)[1]))

		cl.ch = rbind(pm.cl, ponly, monly)
		cl.ch$chr = c
		cl.ch$size = cl.ch$e - cl.ch$s + 1
		cl.ch$ram = cl.ch$ram1 + cl.ch$ram2
		cl.ch$rep.repeat = unlist(lapply(1:dim(cl.ch)[1], function(i) {
		   return(paste(unique(c(strsplit(as.character(cl.ch$rep.repeat1[i]), ",")[[1]], strsplit(as.character(cl.ch$rep.repeat2[i]), ",")[[1]])), collapse=","))
		}))
		cl.ch$family = unlist(lapply(1:dim(cl.ch)[1], function(i) {
		  if (cl.ch$family1[i] == "") return(cl.ch$family2[i]) else if
		     (cl.ch$family2[i] == "") return(cl.ch$family1[i]) else
		    return(cl.ch$family1[i]) # otherwise family1 is equal to family2  
		}))
		cl.ch$class = unlist(lapply(1:dim(cl.ch)[1], function(i) {
		  if (cl.ch$class1[i] == "") return(cl.ch$class2[i]) else if
		     (cl.ch$class2[i] == "") return(cl.ch$class1[i]) else
		    return(cl.ch$class1[i]) # otherwise family1 is equal to family2 
		}))

		# write the raw clusters
		cl.ch = cl.ch[, colh$cl.raw]
		write.table(cl.ch, cl.raw.file, sep="\t", row.names=F, col.names=F, quote=F, append=T)
		write.msg(paste("done writing", dim(cl.ch)[1], "raw clusters for", c))

		colnames(cl.ch)[9:14] = c("pram", "nram", "pram_start", "pram_end", "nram_start", "nram_end")

		# remove clusters whose boundary size is smaller than 2 * rl; related with update on the count.clipped2
        cl.ch = cl.ch[cl.ch$size >= 2 * rl, ]

		# filter out the clusters with less than minimum rams
		cl.ch = cl.ch[cl.ch$ram >= min.ram, ]
        if (nrow(cl.ch) == 0 ) {
            write.msg(sprintf("no cluster with %d rams", min.ram))
            return(NULL)
        }
        if (!oneside.ram) cl.ch = cl.ch[cl.ch$pram >= 1 & cl.ch$nram >= 1, ]

		# count clipped reads
		if (no.clipped) {
			 clipped.cnt = data.frame(acr=0, pacr=0, nacr=0, cr=0, pcr=0, ncr=0, pbp=-1, nbp=-1, tsd=-9999, arr=0)
        } else {
          write.msg("counting clipped reads ...")
          clipped = count.clipped2(cl.ch, cbam.file, verbose, jittering, cbam.chr, rl=rl)
          clipped.cnt = data.frame(do.call(rbind, lapply(clipped, function(x) unlist(x[["ccnt"]]))))
          write.msg(paste("done:", sum(clipped.cnt$cr), "clipped reads"))
        }
        # for compatibility, update the clipped.cnt column header 
        colnames(clipped.cnt) = c("acr", "pacr", "nacr", "cr", "pcr", "ncr", "pbp", "nbp", "tsd", "acrr")
        cl.ch = cbind(cl.ch, clipped.cnt)

		# calculate scores
		cl.ch$score = get.score(cl.ch,  20, 10)

		# annotate the repeat masker te instances

		# mark clusters near known instances: new labels in ril =>  te_type:class:family
		rtype = sapply(strsplit(names(ril), ":"), function(x) x[[1]])
		rclass = sapply(strsplit(names(ril), ":"), function(x) x[[2]])
		rfamily = sapply(strsplit(names(ril), ":"), function(x) x[[3]])
		# note that the rmasker family for SVA is "OTHER"
		rfamily[grep("OTHER", rfamily)] = "SVA"	
	
		write.msg("marking clusters near known te instances ...")
		ins.margin = is$ins.margin
		if (ins.margin < rmasker.filter.margin) {
				ins.margin = rmasker.filter.margin 
				print(sprintf("the estimated ins.margin is too small:%d.", is$ins.margin))
		}

		cl.ch$oi = "-" # not annotated
		if (annot.oi) {
			write.msg("annotating known repeat instances ...")
 			urepeats = unique(cl.ch$rep.repeat)
 			names(urepeats) = urepeats 
			# prepare the merged ril for each repeat type for chromosomec 
			mril = lapply(urepeats, function(r) return(do.call(rbind, lapply(ril[grep(r, rfamily)], function(x) x[[c]]))))
			cl.ch$oi = annotate.oi.chr(cl.ch, mril, rmasker.filter.margin, verbose)
			write.msg("done annotating.")
		}

		# mark insertions near the reference gap regions
		cl.ch$gap=-1 # not defined
		if (annot.gap & !is.null(ref.annot$gap)) {
			# use pbp/nbp for the clusters with clipped reads, otherwise use ram boundaries
			cl.ch$gap = region.overlap.chr(cl.ch[,c("chr", "s", "e")], ref.annot$gap[[c]], margin=5)
			idx = which(cl.ch$pbp !=-1 & cl.ch$nbp !=-1)
			cl.ch$gap[idx] = region.overlap.chr(data.frame(chr=cl.ch$chr[idx], s=pmin(cl.ch$pbp[idx], cl.ch$nbp[idx]), e=pmax(cl.ch$pbp[idx], cl.ch$nbp[idx])), ref.annot$gap[[c]], margin=5)
		}

        cl.ch$pgene = cl.ch$ngene = "-"
        if (annot.gene) {
            write.msg("annotating genes ...")
            cl.ch = annotate.genes(cl.ch, ref.annot$genes, margin=2)
            write.msg("done annotating.")
        }

		if (exo && rasym != "um") cl.ch$desc = annotate.virus(cl.ch) else cl.ch$desc = "-"
		cl.ch = cl.ch[, colh$cl]
		write.table(cl.ch, cl.file, sep="\t", quote=F, row.names=F, append=T, col.names=F)
		write.msg(paste("done writing", nrow(cl.ch), "clusters for", c))

		if (no.clipped) {
			clipped.detail=NULL
		} else {
			clipped.detail = data.frame(do.call(rbind, lapply(clipped, function(x)
		 	 if (!is.null(x$cdf)) return(x[["cdf"]]))))
			if (nrow(clipped.detail)>0) {
	        	colnames(clipped.detail) = c("chr", "s", "e", "cpos", "cigar", "tname", "ref.seq", "clipped.seq", "clipped.qual", "aligned")
		 	 	clipped.detail = clipped.detail[, colh$clipped]
		 	 	write.table(clipped.detail, clipped.file,  row.names=F, quote=F, sep="\t", append=T, col.names=F)
		 	 	write.msg(paste("done writing", nrow(clipped.detail), "clipped sequence for", c))
			}
		}
		rm(cl.ch, clipped.detail)

		write.msg(paste(c, "elapsed:", (proc.time()-ptm)[["elapsed"]]))

		return(1)
	})
	rm(cl, ram)
	write.msg(paste("all elapsed:", (proc.time()-ptm0)[["elapsed"]]))

	return(1)
}


# the original tea
rid.te <- function(sample, dir, chr=NULL, no.clipped=F, ref="hg18", min.ram=3, jittering=2, cbam.chr=T, verbose=F)
{

	#rilf= "/groups/park/alee/ra/data/rmasker/hg18.allri.merged.repeat.clusters.txt.RData"
	rilf = sprintf("%s/lib/%s.allri.merged.repeat.clusters.txt.RData", tea.base, ref)

	bam.dir = paste(dir, "/", sample, "/bam/",  sep="")
	cl.dir = paste(dir, "/", sample, "/cluster/",  sep="")
	if (!file.exists(cl.dir)) { dir.create(cl.dir, recursive=T, mode="0755")}

	if (is.null(chr)) out.prefix = paste(cl.dir, sample, sep="") else
	  out.prefix = paste(cl.dir, sample, ".", chr, sep="")

	ram.file = paste(bam.dir, sample, ".ram.bz2", sep="")
 
	cbam.file = NULL
	if (!no.clipped) cbam.file = paste(bam.dir, sample, ".sorted.softclips.consd.bam", sep="")

	rl.file = paste(bam.dir, sample, ".rl", sep="")
	isize.file = paste(bam.dir, sample, ".isize", sep="")

	write.msg(paste("loading rl and isize from", rl.file, isize.file, "..."))
	rl = load.rl(rl.file)
	is = load.isize(isize.file, rl)
	write.msg(paste("fragment:", is$fr, "(mu:", is$mu, "sd:", is$sd, ") intra.gap:", is$intra.gap, " inter.gap:", is$inter.gap, " ins.margin:", is$ins.margin))

	write.msg(paste("loading rams from", ram.file, "..."))
	ram = load.ram(ram.file, separated=F)
	write.msg("done loading rams.")

	if (!is.null(chr)  && is.null(ram[[chr]]))  { stop(sprintf("there is no ram for %s", chr)) }

	write.msg(paste("loading repeat instances from", rilf, "..."))
	load(rilf)
	write.msg("done loading instances.")

	if (is.null(chr))   chrl = intersect(names(ram), paste("chr", c(1:22, "X", "Y"), sep="")) else
	  chrl = c(chr)
	names(chrl) = chrl

	# write headers
    col.cl.raw = c("chr", "s", "e", "size", "rep.repeat", "family", "class", "ram", "ram1", "ram2", "s1", "e1",
            "s2", "e2", "pos1", "pos2", "rep.repeat1", "rep.repeat2", "repeat.name1", "repeat.name2", "rname1", "rname2")
    col.cl = c("chr", "s", "e", "size", "td", "cp1", "cp2", "rep.repeat", "family", "class", "ram", "ram1", "ram2",
            "cr", "cr1", "cr2", "acr", "acr1", "acr2", "arr", "s1", "e1", "s2", "e2", "score", "oi")
    col.clipped = c("chr", "s", "e", "cs", "ce", "cpos", "aligned", "cigar", "tname", "ref.seq", "clipped.seq", "clipped.qual")

	if (is.null(chr) || chr == "chr1") {
	  write(paste(col.cl.raw, collapse="\t"), paste(out.prefix, ".cluster.raw", sep=""))
	  write(paste(col.cl, collapse="\t"), paste(out.prefix, ".cluster", sep=""))
	  write(paste(col.clipped, collapse="\t"), paste(out.prefix, ".clipped", sep=""))
	}

	ptm0 = proc.time()
	cl = lapply(chrl, function(c) {
	  write.msg(paste("processing", c))
	  ptm = proc.time()
	  if (dim(ram[[c]]$p)[1] == 0 | dim(ram[[c]]$m)[1] == 0) return (NULL)

	  write.msg(paste("extracting positive strand ram (pram) clusters ..."))
	  ram[[c]]$p$repeatname = gsub("/", "_", ram[[c]]$p$repeatname) # replace ALR/ALPHA to ALR_ALPHA
	  p.cl = get.cluster(ram[[c]]$p$pos, 1, ram[[c]]$p$repeatname, ram[[c]]$p$rname, is$intra.gap, merge.family, verbose)
	  write.msg(paste("done: positive strand clusters:", dim(p.cl)[1]))

	  write.msg(paste("extracting negative strand  ram (nram) clusters ..."))
	  ram[[c]]$m$repeatname = gsub("/", "_", ram[[c]]$m$repeatname) # replace ALR/ALPHA to ALR_ALPHA
	  m.cl = get.cluster(ram[[c]]$m$pos, -1, ram[[c]]$m$repeatname, ram[[c]]$m$rname, is$intra.gap, merge.family, verbose)
	  write.msg(paste("done: negative strand clusters:", dim(m.cl)[1]))

	  write.msg(paste("pairing clusters ..."))
	  pov = pair.cluster(p.cl, m.cl, is$inter.gap)
	  pm.cl = pov$pcl
	  write.msg(paste("done: paired clusters:", dim(pm.cl)[1]))

	if (dim(pm.cl)[1] >0) { 
	    # define ram cluster boundaries and remove clusters whose boundary size is smaller than 2 * rl + 10 
	    pm.cl$s = pm.cl$s1
	    pm.cl$e = pm.cl$e2+rl
	    pm.cl = pm.cl[pm.cl$e - pm.cl$s + 1 >= 2 * rl + 10, ]
	
	    colnames(p.cl) = paste(colnames(p.cl), '1', sep="")
	    colnames(m.cl) = paste(colnames(m.cl), '2', sep="")
	
	    empty.len = length(setdiff(1:dim(p.cl)[1], pov$ov[,1]))
	    empty.num = rep(0, empty.len)
	    empty.chr = rep("", empty.len)
	    empty.df  = data.frame(empty.num, empty.num, empty.num, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr)
	    colnames(empty.df) = colnames(m.cl)
	
	    ponly = cbind(p.cl[setdiff(1:dim(p.cl)[1], pov$ov[,1]),], empty.df)
	    ponly$s = ponly$s1
	    ponly$e = ponly$e1 + is$fr
	    write.msg(paste("positive ram only clusters:", dim(ponly)[1]))
	
	    empty.len = length(setdiff(1:dim(m.cl)[1], pov$ov[,2]))
	    empty.num = rep(0,empty.len)
	    empty.chr = rep("", empty.len)
	    empty.df  = data.frame(empty.num, empty.num, empty.num, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr, empty.chr)
	    colnames(empty.df) = colnames(p.cl)
	    monly = cbind(empty.df, m.cl[setdiff(1:dim(m.cl)[1], pov$ov[,2]),])
	    monly$s = monly$s2 - is$fr
	    monly$e = monly$e2
	    write.msg(paste("negative ram only clusters:", dim(monly)[1]))
	
	    cl.ch = rbind(pm.cl, ponly, monly)
	    cl.ch$chr = c
	    cl.ch$size = cl.ch$e - cl.ch$s + 1
	    cl.ch$ram = cl.ch$ram1 + cl.ch$ram2
	    cl.ch$rep.repeat = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	       return(paste(unique(c(strsplit(as.character(cl.ch$rep.repeat1[i]), ",")[[1]], strsplit(as.character(cl.ch$rep.repeat2[i]), ",")[[1]])), collapse=","))
	    }))
	    cl.ch$family = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	      if (cl.ch$family1[i] == "") return(cl.ch$family2[i]) else if
	         (cl.ch$family2[i] == "") return(cl.ch$family1[i]) else
	        return(cl.ch$family1[i]) # otherwise family1 is equal to family2  
	    }))
	    cl.ch$class = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	      if (cl.ch$class1[i] == "") return(cl.ch$class2[i]) else if
	         (cl.ch$class2[i] == "") return(cl.ch$class1[i]) else
	        return(cl.ch$class1[i]) # otherwise family1 is equal to family2 
	    }))
	
	    cl.ch = cl.ch[, col.cl.raw]
	    write.table(cl.ch, paste(out.prefix, ".cluster.raw", sep=""), sep="\t", row.names=F, col.names=F, quote=F, append=T)
	    write.msg(paste("done writing", dim(cl.ch)[1], "raw clusters for", c))
	
	  # filter out the clusters with less than minimum rams
	    cl.ch = cl.ch[cl.ch$ram >= min.ram & cl.ch$ram1 >= 1  & cl.ch$ram2 >= 1, ]
	
	    if (no.clipped) {
	      clipped.cnt = data.frame(td=-9999, cp1=-1,cp2=-1, cr=0, cr1=0, cr2=0, acr=0, acr1=0, acr2=0, arr=0)
	    } else {
	      write.msg("counting clipped reads ...")
	      clipped = count.clipped2(cl.ch, cbam.file, verbose, jittering, cbam.drop.chr, rl=rl)
	      clipped.cnt = data.frame(do.call(rbind, lapply(clipped, function(x) unlist(x[["ccnt"]]))))
	      write.msg(paste("done:", sum(clipped.cnt$cr), "clipped reads"))
	    }
			
	   cl.ch = cbind(cl.ch, clipped.cnt)
	
	    # calculate scores
	    cl.ch$score = get.score(cl.ch[, c("ram1", "ram2", "acr1", "acr2")], 20, 10)
	
	    # mark clusters near known instances
	    write.msg("marking clusters near known repeat instances ...")
	    cl.ch$oi = unlist(lapply(1:dim(cl.ch)[1], function(i) {
	      if (verbose & i%%100 ==0) print(paste("processing", i))
	      repeats = strsplit(cl.ch$rep.repeat[i], ",")[[1]]
	      oi = min(unlist(lapply(repeats, function(r) {
	        oi = ifelse(countOverlaps(IRanges(cl.ch$s[i], cl.ch$e[i]), IRanges(ril[[r]][[c]]$s, ril[[r]][[c]]$e) + is$ins.margin)==0, 1, 0)
	      })))
	      return(oi)
	    }))
		write.msg("done annotation.")
	
	    cl.ch = cl.ch[, col.cl]
	    write.table(cl.ch, paste(out.prefix, ".cluster", sep=""), sep="\t", quote=F, row.names=F, append=T, col.names=F)
	    write.msg(paste("done writing", dim(cl.ch)[1], "clusters for", c))
	
	    if (no.clipped) {
	      clipped.detail=NULL
	    } else {
	      clipped.detail = data.frame(do.call(rbind, lapply(clipped, function(x)
	       if (!is.null(x$cdf)) return(x[["cdf"]]))))
	      if (dim(clipped.detail)[1] >0) {
	        colnames(clipped.detail) = c("chr", "s", "e", "cs", "ce", "cpos", "cigar", "tname", "ref.seq", "clipped.seq", "clipped.qual", "aligned")
	        clipped.detail = clipped.detail[, col.clipped]
	        write.table(clipped.detail, paste(out.prefix, ".clipped", sep=""),  row.names=F, quote=F, sep="\t", append=T, col.names=F)
	        write.msg(paste("done writing", dim(clipped.detail)[1], "clipped sequence for", c))
	      }
	    }
	    write.msg(paste(c, "elapsed:", (proc.time()-ptm)[["elapsed"]]))
	    return(cl.ch)
		} else {
			return(NULL)
		}
	})
	write.msg(paste("all elapsed:", (proc.time()-ptm0)[["elapsed"]]))
	return(1)
}

annot.conf <- function(df, oneside.ram, min.ram, min.acr, min.acrr, min.tsd, max.tsd, min.score, ram.cutoff, min.out.conf, no.oi)
{
    if (no.oi) df = df[df$oi==1 | df$oi == "-", ]

	# for compatibility with the older format
	if ("td" %in% colnames(df)) df = rename.column(df)

    df = df[df$ram>=min.ram,]
    if (!oneside.ram) df = df[df$pram>0 & df$nram>0,]
    df$conf = 1

	if (oneside.ram) { 
    	df$conf[df$acr >= min.acr] = 2
	} else {
    	df$conf[df$acr >= min.acr & df$pacr >= floor(min.acr/2) & df$nacr >=floor(min.acr/2)] = 2
	}
    df$conf[df$conf==2 & df$acrr >= min.acrr] = 3

    df$conf[df$conf==3 & df$tsd >=min.tsd & df$tsd <= max.tsd] = 4
    df$conf[df$conf==4 & df$ram>=ram.cutoff] = 5

    df = df[df$conf >=min.out.conf,]
	
	return(df)
}

rename.column <- function(df)
{
	colnames(df)[grep("td", colnames(df))] =  "tsd"
	colnames(df)[grep("ram1", colnames(df))] =  "pram"
	colnames(df)[grep("ram2", colnames(df))] =  "nram"
	colnames(df)[grep("acr1", colnames(df))] =  "pacr"
	colnames(df)[grep("acr2", colnames(df))] =  "nacr"
	colnames(df)[grep("arr", colnames(df))] =  "acrr"

	return(df)
}

# mainly used for primate genomes
call.germline <- function(dir, sample, ref="hg18", rasym="ra", min.ram=3, oneside.ram=F, min.acr=2, min.acrr=0.4, min.tsd=-20, max.tsd=50, no.oi=T, min.out.conf=5, ram.cutoff=6, verbose=F, mark.exo=F, parallel=F, contig=F, seqid=F)
{
	
	ref.annot = load.ref.info(ref)
	chrl = ref.annot$chrl
	# load clusters
	if (is.null(rasym)) cldir = sprintf("%s/%s/cluster", dir, sample) else
		cldir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)

	#write.msg(paste("done generatingYOYO", chrl))
	
	#cl = load.cl.rfile(cldir, sample, chrl, verbose)
	cl = load.cl.rfile(cldir, sample, verbose)
	
	#write.msg(paste("done gener",cl))
    df = data.frame(cl)
    #df = data.frame(do.call(rbind, cl))
    #df=cl
    rownames(df) = NULL

	df = annot.conf(df, oneside.ram, min.ram, min.acr,  min.acrr, min.tsd, max.tsd, min.score, ram.cutoff, min.out.conf, no.oi)
	print(nrow(df))

	df = cbind(sample=sample, df)

	if (is.null(rasym)) outdir = sprintf("%s/%s/cluster", dir, sample) else
    	outdir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)

	out.file = sprintf("%s/%s.germline", outdir, sample)

    write.table(df, out.file, sep="\t", quote=F, row.names=F)
    write.msg(paste("done writing", out.file, nrow(df), "insertions."))

	if (!is.null(rasym) && rasym == "um" && mark.exo) {
		df.new = count.ram(df, dir, sample, chrl, fr.ram=0.05, max.ram=5, margin=0)
		exo = df.new[df.new$endo.ram <= df.new$ram.cutoff,]
    	exo.out.file = sprintf("%s/%s.germline.exo", outdir, sample)
	    write.table(exo, exo.out.file, sep="\t", quote=F, row.names=F)
   	 	write.msg(paste("done.writing", out.file, nrow(exo), "insertions."))
	}

	if (contig) {
		exo.f = ""
		parallel.f = ""
		if (rasym != "ra") exo.f = "-x"
		if (parallel) parallel.f = "-p"
		cmd = sprintf("%s/scripts/tea.pl contig -o %s -d %s %s %s -R %s %s", tea.base, outdir, dir, exo.f, parallel, rasym, out.file)
		print(cmd)
		system(cmd)
	}

	#if (seqid && rasym == "um") {
	#	cmd = sprintf("%s/scripts/seqid.v2.sh %s.contig", tea.base, out.file)
	#	print(cmd)
	#	system(cmd)
	#}

    return(1)
}

count.ram <- function(df, dir, sample, chrl, fr.ram=0.05, max.ram=5, margin=0, verbose=F)
{
	# load ram counts
	cldir = sprintf("%s/%s/cluster_ram", dir, sample)
	cldir2 = sprintf("%s/%s/cluster", dir, sample) # for compatibility

	if (!file.exists(cldir) && !file.exists(cldir2)) {
		sprintf("no ram cluster dir: %s, %s", cldir, cldir2)
		return(NULL)
	}

	if (!file.exists(cldir) && file.exists(cldir2)) cldir = cldir2


	ram.raw.rfile = sprintf("%s/%s.cluster.raw.RData", cldir, sample)

	if (file.exists(ram.raw.rfile)) ram = readRDS(ram.raw.rfile)  else
	ram = make.ram.raw.rfile(dir, sample, ram.raw.rfile, chrl, verbose, rasym = "ra")

	chrl = intersect(chrl, unique(df$chr))
	df.new = data.frame(do.call(rbind, lapply(chrl, function(chr) {
		df.ch = df[df$chr==chr,]
		if (nrow(df.ch)==0) return(NULL)

		ram.ch = do.call(rbind, ram[[chr]])
		rownames(ram.ch) = NULL
		if (nrow(ram.ch)==0) return(NULL)
		r1 = IRanges(df.ch$s-margin, df.ch$e+margin)
		r2 = IRanges(as.numeric(ram.ch$s), as.numeric(ram.ch$e))
		# to use the newer version of IRanges
    	#ovm = findOverlaps(r1, r2)@matchMatrix
    	ovm = as.matrix(findOverlaps(r1, r2))

		df.ch$endo.ram = 0

		if (nrow(ovm) > 0) {
        	cram = do.call(rbind, lapply(unique(ovm[,1]), function(i) {
          	return(c(i, max(ram.ch$ram[ovm[ovm[,1] == i,2]])))
			}))
  			df.ch$endo.ram[cram[,1]] = cram[,2]
		}

		return(df.ch)
	})))
	rownames(df.new) = NULL
	df.new$ram.cutoff = pmin(round(df.new$ram * fr.ram), max.ram)	
	
	return(df.new)
}

count.cram <- function(cl, cl.rawl, margin=100, verbose=F)
{
    cl.new = data.frame(do.call(rbind, lapply(unique(cl$family), function(f) {
      if (verbose) print(paste("processing", f))
      cl.f = cl[cl$family==f,]
      cl.f$cram = 0
      if (!is.null(cl.rawl[[f]])) {
        control = cl.rawl[[f]]
		# to use the newer version of IRanges
        #ovm = findOverlaps(IRanges(cl.f$s-margin, cl.f$e+margin), IRanges(as.numeric(control$s), as.numeric(control$e)))@matchMatrix
        ovm = as.matrix(findOverlaps(IRanges(cl.f$s-margin, cl.f$e+margin), IRanges(as.numeric(control$s), as.numeric(control$e))))
        cram = do.call(rbind, lapply(unique(ovm[,1]), function(i) {
          return(c(i, max(control$ram[ovm[ovm[,1] == i,2]])))
        }))
        cl.f$cram[cram[,1]] = cram[,2]
      }
      return(cl.f)
    })))
    cl.new = cl.new[order(as.numeric(rownames(cl.new))),]
    return(cl.new$cram)
}

call.somatic <- function(dir, sample, rasym="ra", matched.control=NULL, nonmatched.controls=NULL, ref="hg18", oneside.ram=F, min.ram=3, min.acr=2, min.acrr=0.5, min.tsd=-15, max.tsd=30, min.score=0.6, verbose=F, gene.annot=T, ram.cutoff=6, matched.cram=1, nonmatched.cram=2, matched.cacr=1, nonmatched.cacr=2, min.out.conf=5, no.oi=T, mark.exo=T, parallel=F, contig=T, seqid=T)
{
	ref.annot = load.ref.info(ref)
	chrl = ref.annot$chrl

	controls = matched.control
	if (!is.null(nonmatched.controls)) controls = c(matched.control, nonmatched.controls)
	names(controls) = controls

	cldir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)
	cl = load.cl.rfile(cldir, sample, chrl, verbose)

	if (rasym == "va") cl = split(cl, cl$chr)
	cl = lapply(cl, function(df) annot.conf(df=df, oneside.ram, min.ram, min.acrr, min.acrr, min.tsd, max.tsd, min.score, ram.cutoff, min.out.conf, no.oi))

	cram  = load.cram.rfile(dir, sample, rasym, cl, controls, chrl, verbose)
	ccr = load.ccr.rfile(dir, sample, rasym, cl, controls, chrl, verbose)

	cl.df = data.frame(do.call(rbind, cl))
	cram.df = data.frame(do.call(rbind, cram))
	ccr.df = data.frame(do.call(rbind, ccr))
	rownames(cl.df) = rownames(cram.df) = rownames(ccr.df) = NULL
	#write.msg(paste(cl.df))
	#write.msg(paste(cram.df))

	cl.df = count.controls(cl.df, cram.df, ccr.df, matched.control, nonmatched.controls)
	cl.df = cbind(sample=sample, cl.df)

	outdir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)
	raw.out.file = sprintf("%s/%s/cluster_%sm/%s.somatic.raw", dir, sample, rasym, sample)
	m.out.file = sprintf("%s/%s/cluster_%sm/%s.somatic.matched", dir, sample, rasym, sample)
	out.file = sprintf("%s/%s.somatic", outdir, sample)
#    exo.out.file = sprintf("%s/%s/cluster_%sm/%s.somatic.exo", dir, sample, rasym, sample)

	write.table(cl.df, raw.out.file, sep="\t", quote=F, row.names=F)
    write.msg(paste("done.writing", raw.out.file, nrow(cl.df), "insertions."))

	cl.df= cl.df[cl.df$m.cram <= matched.cram & cl.df$m.cacr <= matched.cacr,]
	write.table(cl.df, m.out.file, sep="\t", quote=F, row.names=F)
    write.msg(paste("done.writing", m.out.file, nrow(cl.df), "insertions."))

	cl.df= cl.df[cl.df$m.cram <= matched.cram & cl.df$m.cacr <= matched.cacr & cl.df$max.nm.cram <= nonmatched.cram & cl.df$max.nm.cacr <= nonmatched.cacr,]
	write.table(cl.df, out.file, sep="\t", quote=F, row.names=F)
    write.msg(paste("done.writing", out.file, nrow(cl.df), "insertions."))

#	if (rasym == "um" && mark.exo) {
#    	df.new = count.ram(cl.df, dir, sample, chrl, fr.ram=0.2, max.ram=5, margin=0)
#        exo = df.new[df.new$endo.ram <= df.new$ram.cutoff,]
#        write.table(exo, exo.out.file, sep="\t", quote=F, row.names=F)
#        write.msg(paste("done.writing", exo.out.file, nrow(exo), "insertions."))
#	}

	if (contig) {
		exo.f = ""
		parallel = ""
		if (rasym != "ra" && rasym != "youngte") exo.f = "-x"
		#if (parallel) parallel.f = "-p"

		#if (file.exists(sprintf("%s/%s.cluster.bak", outdir, sample))) {
			cmd = sprintf("mv %s/%s.cluster.bak %s/%s.cluster", outdir, sample, outdir, sample)
			system(cmd)
		#}
		#if (file.exists(sprintf("%s/%s.cluster.raw.bak", outdir, sample))) {
			cmd = sprintf("mv %s/%s.cluster.raw.bak %s/%s.cluster.raw", outdir, sample, outdir, sample)
			system(cmd)
		#}
		#if (file.exists(sprintf("%s/%s.clipped.bak ", outdir, sample))) {
			cmd = sprintf("mv %s/%s.clipped.bak %s/%s.clipped", outdir, sample, outdir, sample)
			system(cmd)
		#}

		parallel.f = 1
		cmd = sprintf("%s/scripts/tea.pl contig -o %s -d %s %s -R %s %s", tea.base, outdir, dir, exo.f, rasym, out.file)
		#write.msg(paste(cmd))
		#print(cmd)
		system(cmd)
	}

	if (seqid && rasym == "um") {
		cmd = sprintf("%s/scripts/seqid.v2.sh %s.contig", tea.base, out.file)
		#print(cmd)
		system(cmd)
	}

	return(1)
}

count.controls <- function(cl.df, cram.df, ccr.df, matched.control, nonmatched.controls)
{
	# to handle the automatic colume name conversion of a data frame
	# add "X" in the column names if the sample name starts with a number
	
	matched.control[grep("^[0-9].*", matched.control)] = paste("X", matched.control[grep("^[0-9].*", matched.control)], sep="")
	matched.control = gsub("-", ".", matched.control)
	
	m.cram.idx = match(matched.control, colnames(cram.df))
	#write.msg(paste(m.cram.idx))
	m.ccr.idx = match(matched.control, colnames(ccr.df))
	
	cl.df$m.cram = cram.df[, m.cram.idx]
	cl.df$m.cacr = ccr.df[, m.ccr.idx]

	if (is.null(nonmatched.controls)) {
		cl.df$max.nm.cram = 0  
    	cl.df$max.nm.cacr = 0 

    	cl.df$nm.cram.cnt = "-"
    	cl.df$nm.cacr.cnt = "-"
	} else {
		nonmatched.controls[grep("^[0-9].*", nonmatched.controls)] = paste("X", nonmatched.controls[grep("^[0-9].*", nonmatched.controls)], sep="")
		nonmatched.controls = gsub("-", ".", nonmatched.controls)
		nm.cram.idx = match(nonmatched.controls, colnames(cram.df)) # sometimes there are missing controls in cram.df
		
		#write.msg(paste(colnames(cram.df)))
		#write.msg(paste(nm.cram.idx))
		nm.cram.cnames = colnames(cram.df)[nm.cram.idx]

		nm.ccr.idx = match(nonmatched.controls, colnames(ccr.df))
		nm.ccr.cnames = colnames(ccr.df)[nm.ccr.idx]
		
		#write.msg(paste(cram.df[, nm.cram.idx]))
		cl.df$max.nm.cram = apply(cram.df[, nm.cram.idx], 1, max)
		
		cl.df$max.nm.cacr = apply(ccr.df[, nm.ccr.idx], 1, max)

		cl.df$nm.cram.cnt = apply(cram.df[, nm.cram.idx], 1, function(x)
            paste(paste(nm.cram.cnames[x!=0], x[x!=0], sep=":"), collapse=","))

		cl.df$nm.cacr.cnt = apply(ccr.df[, nm.ccr.idx], 1, function(x)
       		paste(paste(nm.ccr.cnames[x!=0], x[x!=0], sep=":"), collapse=","))
	}

    return(cl.df)
}

write.somatic.matchc.df2 <- function(cl.df, cram.df, ccr.df, controls, sample, max.cram, max.ccr, fprefix, verbose=T)
{
	cram.idx = match(controls, colnames(cram.df))
	cram.cnames = colnames(cram.df)[cram.idx]

	ccr.idx = match(controls, colnames(ccr.df))
	ccr.cnames = colnames(ccr.df)[ccr.idx]

	if (grepl("cancer", sample))  {
	  cram.cidx = match(sub("cancer", "normal", sample), colnames(cram.df))
	  ccr.cidx = match(sub("cancer", "normal", sample), colnames(ccr.df))
	} else if (grepl("normal", sample))  { # for extracting normal-specific insertions
	  cram.cidx = match(sub("normal", "cancer", sample), colnames(cram.df))
	  ccr.cidx = match(sub("normal", "cancer", sample), colnames(ccr.df))
	}

	cram.cidx=cram.idx=1; ccr.cidx=ccr.idx=1
	idx = which(cram.df[,cram.cidx] <= max.cram & ccr.df[,ccr.cidx] <= max.ccr)
	cram.cnt = apply(cram.df[idx,cram.idx,drop=FALSE], 1, function(x)
	  paste(paste(cram.cnames[x!=0], x[x!=0], sep=":"), collapse=","))
	ccr.cnt = apply(ccr.df[idx,ccr.idx,drop=FALSE], 1, function(x)
	  paste(paste(ccr.cnames[x!=0], x[x!=0], sep=":"), collapse=","))

	df = cbind(cl.df[idx,,drop=FALSE], cram.cnt, ccr.cnt)
	#df = col.format(df, sub("_(cancer|normal)", "", sample), type)
	df = col.format(df, sub("_(cancer|normal)", "", sample))

	fname = paste(fprefix, max.cram, "cram.", max.ccr, "ccr", sep="")
	write.table(df, file=fname, sep="\t", quote=F, row.names=F)
	write.msg(paste("done.writing", fname))

	return(list(df=df, fname=fname))
}


load.cl.rfile <- function(cldir, sample, chrl, verbose)
{
	cl.rfile = sprintf("%s/%s.cl.RData", cldir, sample)
	if (file.exists(cl.rfile)) {
	  cl = readRDS(cl.rfile) # cl
	  write.msg(paste("done loading", cl.rfile))
	} else {
	  write.msg(paste("generating", cl.rfile))
	  cl = make.cl.rfile2(cldir, sample, cl.rfile, chrl, verbose)
	  write.msg(paste("done generating", cl.rfile))
	}
	return(cl)
}

load.cram.rfile <- function(dir, sample, rasym, cl, controls, chrl, verbose)
{
	cram.rfile = sprintf("%s/%s/cluster_%sm/%s.cram.RData", dir, sample, rasym, sample)
	if (file.exists(cram.rfile)) {
		cram = readRDS(cram.rfile)
		write.msg(paste("done loading", cram.rfile))
	} else {
		write.msg(paste("generating", cram.rfile))
		cram = make.cram.rfile(dir, sample, cl, controls, cram.rfile, chrl, verbose, rasym)
	 	write.msg(paste("done generating", cram.rfile))
	}
	return(cram)
}

load.ccr.rfile <- function(dir, sample, rasym, cl, controls, chrl, verbose)
{
	ccr.rfile = sprintf("%s/%s/cluster_%sm/%s.ccr.RData", dir, sample, rasym, sample)
	if (file.exists(ccr.rfile)) {
		ccr = readRDS(ccr.rfile) 
		write.msg(paste("done loading", ccr.rfile))
	} else {
		write.msg(paste("generating", ccr.rfile))
		ccr = make.ccr.rfile (dir, sample, cl, controls, ccr.rfile, chrl, verbose)
		write.msg(paste("done generating", ccr.rfile))
	}
	return(ccr)
}

make.cram.rfile <- function(dir, sample, cl, controls, fn, chrl, verbose=F, rasym=NULL)
{
	cram = lapply(chrl, function(chr) {
	  if (verbose) print(paste("processing", chr))
	  return(data.frame(do.call(cbind, lapply(controls, function(control) {
	    if (verbose) print(paste("processing", control))
	      get.cram(dir, control, cl[[chr]], chr, 100, chrl, verbose, rasym)
	  }))))
	})
	saveRDS(cram, file=fn)
	write.msg(paste("done writing the control ram file",fn))
	return(cram)
}

get.cram <- function(dir, control, cl.ch, chr, margin=100, chrl, verbose, rasym=NULL)
{
	# load raw ram counts in a control sample
	if (!is.null(rasym)) {
		ram.raw.rfile = sprintf("%s/%s/cluster_%sm/%s.cluster.raw.RData", dir, control, rasym, control)
	} else {
		ram.raw.rfile = sprintf("%s/%s/cluster/%s.cluster.raw.RData", dir, control, control)
	}

	if (file.exists(ram.raw.rfile)) {
		sprintf("loading %s", ram.raw.rfile)
		raw = readRDS(ram.raw.rfile)
	} else {
		sprintf("making %s", ram.raw.rfile)
		raw = make.ram.raw.rfile(dir, control, ram.raw.rfile, chrl, verbose, rasym)
	}

	cram.chr = count.cram(cl.ch, raw[[chr]], margin=margin, verbose=verbose)

	return(cram.chr)
}

make.ram.raw.rfile <- function(dir, sample, fn, chrl, verbose, rasym)
{
	cols = c("chr", "s", "e", "rep.repeat", "family", "class", "ram", "ram1", "ram2")
	raw = lapply(chrl, function(chr) {
		if (verbose) print(paste("processing", chr))
		if (!is.null(rasym)) {
			raw.cl.file = sprintf("%s/%s/cluster_%sm/%s.%s.cluster.raw", dir, sample, rasym, sample, chr)
		} else {
			raw.cl.file = sprintf("%s/%s/cluster/%s.%s.cluster.raw", dir, sample, sample, chr)
		}
    	cl = read.delim(raw.cl.file, sep="\t", as.is=T)
	  	cl = cl[, cols]
		familyl = unique(cl$family); names(familyl) = familyl
		cl.raw = lapply(familyl, function(f) return(cl[cl$family==f,]))

		return(cl.raw)
	})
	saveRDS(raw, file=fn)
	write.msg(paste("done writing rawl for", sample, ":", fn))
	return(raw)
}

count.cram <- function(cl, cl.rawl, margin=100, verbose=F)
{
	cl.new = data.frame(do.call(rbind, lapply(unique(cl$family), function(f) {
	  if (verbose) print(paste("processing", f))
	  cl.f = cl[cl$family==f,]
	  cl.f$cram = 0
	  if (!is.null(cl.rawl[[f]])) {
	    control = cl.rawl[[f]]
	    ovm = as.matrix(findOverlaps(IRanges(cl.f$s-margin, cl.f$e+margin), IRanges(as.numeric(control$s), as.numeric(control$e))))
	    cram = do.call(rbind, lapply(unique(ovm[,1]), function(i) {
	      return(c(i, max(control$ram[ovm[ovm[,1] == i,2]])))
	    }))
	    cl.f$cram[cram[,1]] = cram[,2]
	  }
	  return(cl.f)
	})))
	cl.new = cl.new[order(as.numeric(rownames(cl.new))),]
	return(cl.new$cram)
}

make.ccr.rfile <- function(dir, sample, cl, controls, fn, chrl, verbose=F)
{
	ccr.t = lapply(controls, function(c) {
	  get.ccr(dir, sample, c, cl, margin=2, chrl, verbose)
	})
	ccr = lapply(chrl, function(chr) {
	  data.frame(do.call(cbind, lapply(controls, function(c) {
	  return(ccr.t[[c]][[chr]])
	  })))
	})
	saveRDS(ccr, file=fn)
	write.msg(paste("done writing the control clipped reads file for sample", sample, ":", fn))
	return(ccr)
}

get.ccr <- function(dir, sample, control, cl, margin, chrl, verbose=F)
{
	cpos.file1 = sprintf("%s/%s/bam/%s.sorted.softclips.consd.cpos.bz2", dir, control, control)
	cpos.file2 = sprintf("%s/%s/bam/%s.softclips.consd.cpos.bz2", dir, control, control)
	if (file.exists(cpos.file1)) {
		cpos.file = cpos.file1 
	} else if (file.exists(cpos.file2)) {
		cpos.file  = cpos.file2
	} else {
		stop(sprintf("no file for clipping positions for sample %s", control))
	}
	cpos.rfile = sub(".bz2", ".RData", cpos.file)

	if (file.exists(cpos.rfile)) {
		write.msg(paste("loading", cpos.rfile))
        cpos = readRDS(cpos.rfile)
	} else {
		write.msg(paste("reading", cpos.file))
		x = scan(cpos.file, what=list(chr='a', pos=1))
		cpos.chrl = unique(as.character(x$chr))
		if (length(grep("chr", cpos.chrl)) == 0) {	
			cpos.chrl = paste("chr", cpos.chrl, sep="")
			drop.chr = T
		} else {
			drop.chr = F
		}
		chrl = intersect(chrl, cpos.chrl)
		names(chrl) = chrl
		cpos = lapply(chrl, function(c) {
	          if (verbose) print(paste("reading", c))
				if (drop.chr) c = sub("chr", "", c)
	          df = data.frame(chr=x$chr[x$chr==c], pos=x$pos[x$chr==c], stringsAsFactors=F)
	          p = sort(df$pos[df$pos > 0])
	          m = sort(-df$pos[df$pos < 0])
	          return(list(p=p, m=m))
	      })
	      saveRDS(cpos, file=cpos.rfile)
	      write.msg(paste("done generating", cpos.rfile))
	}

	ccrl = lapply(chrl, function(chr) {
		cl.ch = cl[[chr]]
	    cl.ch$pbp[cl.ch$pbp==-1] = -1 - margin
	    cl.ch$nbp[cl.ch$nbp==-1] = -1 - margin
	    ccr = points.within(cpos[[chr]]$p, cl.ch$pbp-margin-1, cl.ch$pbp+margin, return.point.counts=T) +
	  	points.within(cpos[[chr]]$m, cl.ch$nbp-margin-1, cl.ch$nbp+margin, return.point.counts=T)
		return(ccr)
	  })
	return(ccrl)
}

call.somatic.old <- function(dir, sample, controls=NULL, min.pram=1, min.nram=1, min.ram=3, min.acr1=1, min.acr2=1, min.arr=0.4, min.tsd=-15, max.tsd=30, min.score=0.6, chrl=hcrl, contig=F, verbose=F)
{
	names(chrl) <- chrl

	for (i in ls())
	  print(paste(i, get(i), sep="="))

	# controls: e.g. controls=c("pr", "ov")
	# using 26 genomes (16 gbm, 8 ov and 2 hapmap) for the non-matched control filtering
	if (is.null(controls)) {
	  #XXX: change me back, uncomment line below
	  #controls = sub("cancer", "normal", sample) # put the matched control in the first 
	  controls = c(controls, dir(dir, "(cr|gbm|ov).*_normal|na18508)"))
	  controls = c(controls, dir(dir, "na1850(6|8)"))
	  # Joe: if there's no occurrence of ov1103_normal in
	  # 'controls', then the grep() returns an empty vector.  Then
	  # the assignment chooses _none_ of the controls.
	  if (sample != "ov1103_cancer") {
	    controls <- controls[grep("ov1103_normal", controls, invert=TRUE)]

	  }
	  #controls = controls[-grep("cra00r_normal|ov0982_normal|pr0508_normal", controls)]
	} else {
	  controls = unlist(sapply(controls, function(t) dir(dir, paste(t, ".*_normal", sep=""))))
	  controls = c(controls, dir(dir, "na1850(7|8)"))
	  # exclude problamatic normal samples
	  controls = controls[grep("cra00r_normal|ov0982_normal|pr0508_normal", controls, invert=TRUE)]
	}
	write.msg("found controls:")
	print(controls)
	names(controls) = controls

	cl.rfile = paste(dir, "/", sample, "/cluster/", sample, ".cl.te.new.RData", sep="")
	cram.rfile = paste(dir, "/", sample, "/cluster/", sample, ".cram.new.RData", sep="")
	ccr.rfile = paste(dir, "/", sample, "/cluster/", sample, ".ccr.new.RData", sep="")

	if (file.exists(cl.rfile)) {
	  load(cl.rfile) # cl
	  write.msg(paste("done loading", cl.rfile))
	} else {
	  write.msg(paste("generating", cl.rfile))
	  cl = make.cl.rfile(dir, sample, cl.rfile, chrl, verbose)
	  write.msg(paste("done generating", cl.rfile))
	}

	if (file.exists(cram.rfile)) {
	  load(cram.rfile)
	  write.msg(paste("done loading", cram.rfile))
	} else {
	  write.msg(paste("generating", cram.rfile))
	    cram = make.cram.rfile(dir, cl, controls, cram.rfile, chrl, verbose)
	  write.msg(paste("done generating", cram.rfile))
	}
	 
	if (file.exists(ccr.rfile)) {
	  load(ccr.rfile)
	  write.msg(paste("done loading", ccr.rfile))
	} else {
	  write.msg(paste("generating", ccr.rfile))
	    ccr = make.ccr.rfile(dir, sample, cl, controls, ccr.rfile, chrl=chrl, verbose=verbose)
	  write.msg(paste("done generating", ccr.rfile))
	}
 
	df = data.frame(do.call(rbind, lapply(chrl, function(chr) {
	  return(cl[[chr]])
	})))
	rownames(df) = NULL

	df$conf=0
	df$conf[df$pram>=1 & df$nram>=1 & df$ram>=3] = 1 # score 0.3
	df$conf[df$conf==1 & df$ram>=6] = 2 # score 0.3
	df$conf[df$conf==2 & df$acr1>=min.acr1 & df$acr2>=min.acr2] = 3 # score 0.5
	df$conf[df$conf==3 & df$tsd >=min.tsd & df$tsd <= max.tsd] = 4
	df$conf[df$conf==4 & df$arr >= min.arr] = 5

	cram.df = data.frame(do.call(rbind, lapply(chrl, function(chr) {
	  return(cram[[chr]])
	})))
	rownames(cram.df) = NULL
	ccr.df = data.frame(do.call(rbind, lapply(chrl, function(chr) {
	  return(ccr[[chr]])
	})))
	rownames(ccr.df) = NULL
	# Joe: drop=FALSE keeps data.frames with only a single column from
	# being coerced to a vector and thus losing the column names.
	cl.df = df[df$conf>=1, , drop=FALSE]
	cram.df = cram.df[df$conf>=1, , drop=FALSE]
	ccr.df = ccr.df[df$conf>=1, , drop=FALSE]

	# apply the matched control filtering 
	out.dir = paste(dir, "/", sample, "/call/", sep="")
	if (!file.exists(out.dir))  dir.create(out.dir, recursive=T, mode="0755")
	fprefix = paste(out.dir, sample, ".somatic.", min.ram, "ram.", min.pram, "pram.", min.nram, "nram.", sep="")

	x = write.somatic.matchc.df2(cl.df, cram.df, ccr.df, controls, sample, max.cram=1, max.ccr=1, fprefix)
	# apply the non-matched control filtering
	x.nmc = write.somatic.nonmatchc.df(x$df, controls, max.cram=2, max.ccr=2, x$fname)
	# annotate conf 
	x.nmc.conf = write.conf(x.nmc$df, min.ram=6, min.acr1=1, min.acr2=1, min.tsd=-15, max.tsd=30, min.arr=0.5, x.nmc$fname)
}

col.format <- function(df, sample, type="somatic")
{
	if (! "gap" %in% colnames(df)) df$gap="-";
	if (type == "somatic") {
	  new.df = data.frame(sample=sample, chr=df$chr, pram_start=df$s, pram_end=df$e1+(df$e[1]-df$e2[1]), nram_start=df$s2, nram_end=df$e, cluster_size=df$size, pbp=df$pbp, nbp=df$nbp, tsd=df$td, rep.repeat=df$rep.repeat, family=df$family, class=df$class, ram=df$ram, pram=df$pram, nram=df$nram, pcr=df$acr1, ncr=df$acr2, acrr=df$arr, score=df$score, oi=df$oi, gap=df$gap, cram=df$cram.cnt, ccr=df$ccr.cnt, stringsAsFactors=F)
	} else if (type == "germline") {
	  new.df = data.frame(sample=sample, chr=df$chr, pram_start=df$s, pram_end=df$e1+(df$e[1]-df$e2[1]), nram_start=df$s2, nram_end=df$e, cluster_size=df$size, pbp=df$pbp, nbp=df$nbp, tsd=df$tsd, rep.repeat=df$rep.repeat, family=df$family, class=df$class, ram=df$ram, pram=df$pram, nram=df$nram, pcr=df$acr1, ncr=df$acr2, acrr=df$arr, score=df$score, oi=df$oi, gap=df$gap, stringsAsFactors=F)
	}
	return(new.df)
}

# using columns before column name formatting
annot.genes <- function(cl, genes, margin=2, verbose=F)
{
	chrl = unique(as.character(cl$chr)); names(chrl) = chrl
	df = data.frame(do.call(rbind, lapply(chrl, function(ch) {
	  if (verbose) print(paste("processing", ch))
	  genes.ch = genes[genes$chr == ch, ]
	  cl.ch = cl[cl$chr==ch,]
	  cl.ch$pgene = "-"
	  #m1 = data.frame(findOverlaps(IRanges(cl.ch$pbp-margin, cl.ch$pbp+margin), IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e)))@matchMatrix)
	  m1 = data.frame(findOverlaps(IRanges(cl.ch$pbp-margin, cl.ch$pbp+margin), IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e)))@matchMatrix)
	  genes = unlist(lapply(unique(m1$query), function(q) {
	    s = m1$subject[m1$query == q]
	    genes = paste(unique(paste(genes.ch[s,]$name, genes.ch[s,]$type, sep="_")), collapse=",")
	    return(genes)
	  }))
	  cl.ch$pgene[unique(m1$query)] = genes

	  cl.ch$ngene = "-"
	  m2 = data.frame(findOverlaps(IRanges(cl.ch$nbp-margin, cl.ch$nbp+margin), IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e)))@matchMatrix)
	  genes = unlist(lapply(unique(m2$query), function(q) {
	    s = m2$subject[m2$query == q]
	    genes = paste(unique(paste(genes.ch[s,]$name, genes.ch[s,]$type, sep="_")), collapse=",")
	    return(genes)
	  }))
	  cl.ch$ngene[unique(m2$query)] = genes
	  return(cl.ch)
	})))
	rownames(df) = NULL
	return(df)
	return(genes)
}

annotate.genes <- function(cl, genes, margin=2, verbose=F)
{
    chrl = unique(as.character(cl$chr)); names(chrl) = chrl
    df = data.frame(do.call(rbind, lapply(chrl, function(ch) {
      if (verbose) print(paste("processing", ch))
      genes.ch = genes[genes$chr == ch, ]
      cl.ch = cl[cl$chr==ch,]
      cl.ch$pgene = "-"
		# to use the newer version of IRanges		
      #m1 = data.frame(findOverlaps(IRanges(cl.ch$pbp, cl.ch$pbp) + margin, IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e)))@matchMatrix)
      #genes = unlist(lapply(unique(m1$query), function(q) {
      #  s = m1$subject[m1$query == q]
      #  genes = paste(unique(paste(genes.ch[s,]$name, genes.ch[s,]$type, sep="_")), collapse=",")
      #  return(genes)
      #}))
      #cl.ch$pgene[unique(m1$query)] = genes
#
#      cl.ch$ngene = "-"
#      m2 = data.frame(findOverlaps(IRanges(cl.ch$nbp, cl.ch$nbp) + margin, IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e)))@matchMatrix)
#      genes = unlist(lapply(unique(m2$query), function(q) {
#        s = m2$subject[m2$query == q]
#        genes = paste(unique(paste(genes.ch[s,]$name, genes.ch[s,]$type, sep="_")), collapse=",")
#        return(genes)
#      }))
#      cl.ch$ngene[unique(m2$query)] = genes
      m1 = as.matrix(findOverlaps(IRanges(cl.ch$pbp, cl.ch$pbp) + margin, IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e))))
      genes = unlist(lapply(unique(m1[,1]), function(q) {
        s = m1[,2][m1[,1] == q]
        genes = paste(unique(paste(genes.ch[s,]$name, genes.ch[s,]$type, sep="_")), collapse=",")
        return(genes)
      }))
      cl.ch$pgene[unique(m1[,1])] = genes

      cl.ch$ngene = "-"
      m2 = as.matrix(findOverlaps(IRanges(cl.ch$nbp, cl.ch$nbp) + margin, IRanges(as.numeric(genes.ch$s), as.numeric(genes.ch$e))))
      genes = unlist(lapply(unique(m2[,1]), function(q) {
        s = m2[,2][m2[,1] == q]
        genes = paste(unique(paste(genes.ch[s,]$name, genes.ch[s,]$type, sep="_")), collapse=",")
        return(genes)
      }))
      cl.ch$ngene[unique(m2[,1])] = genes

      return(cl.ch)
    })))
    rownames(df) = NULL
    return(df)
}

make.cl.rfile2 <- function(cldir, sample, fn, chrl, verbose=F)
{
	cl.file = sprintf("%s/%s.cluster", cldir, sample)

	if (file.exists(cl.file)) {
	  	cl = read.delim(cl.file, sep="\t", as.is=T)
	} else {
		cl = lapply(chrl, function(chr) {
		  if (verbose) print(paste("processing", chr))
			cl.file = sprintf("%s/%s.%s.cluster", cldir, sample, chr)	
			if (!file.exists(cl.file)) {
				print(sprintf("no %s: skipped", cl.file))
				return(NULL) 
			} else {
		  	if (verbose) print(paste("reading", cl.file, sep=""))
		  	x = read.delim(cl.file, sep="\t", as.is=T)
			return(x)
			}
		})
		names(cl) <- chrl
	}
	saveRDS(cl, file=fn)
	return(cl)
}

make.cl.rfile <- function(dir, sample, fn, chrl, verbose=F)
{
    cl = lapply(chrl, function(chr) {
      if (verbose) print(paste("processing", chr))
        cluster.file = sprintf("%s/%s/cluster/%s.%s.cluster", dir, sample, sample, chr)
        if (!file.exists(cluster.file)) {
            print(sprintf("no %s: skipped", cluster.file))
            return(NULL)
        } else {
        	if (verbose) print(paste("reading", cluster.file, sep=""))
        	x = read.delim(cluster.file, sep="\t", as.is=T)
           return(x)
        }
    })
    names(cl) <- chrl
    saveRDS(cl, file=fn)
    return(cl)
}

# seperated: rams are seperated per repeat : three column ram file otherwise four column
load.ram <- function(file, separated=F, rm.dup=F, verbose=F)
{
	if (separated) x = scan(file, what = list(rname='a', chr='a', pos=1)) else
	  x = scan(file, what = list(rname='a', chr='a', pos=1, repeatname='a'))
 
	chrl = unique(as.character(x$chr));
	ram = lapply(chrl, function(c) {
	  if (verbose) print(paste("reading", c))
	  if (separated) df = data.frame(rname = x$rname[x$chr == c], pos = x$pos[x$chr == c], stringsAsFactors=F) else
	    df = data.frame(repeatname = x$repeatname[x$chr==c], pos = x$pos[x$chr == c], rname = x$rname[x$chr == c], stringsAsFactors=F)

	  df.p = df[df$pos > 0,]
	  df.m = df[df$pos < 0,]
	  df.m$pos = -df.m$pos
	  df.p = df.p[order(df.p$pos),]
	  df.m = df.m[order(df.m$pos),]

	  if (rm.dup) {
	    df.p = df.p[!duplicated(df.p$pos),]
	    df.m = df.m[!duplicated(df.m$pos),]
	  }
	  return(list(p=df.p, m=df.m))
	})
	names(ram) = paste("chr", chrl, sep="")

	return(ram)
}

get.cluster2 <- function(sram, strand, gap.cutoff, merge.family=F, verbose=F)
{ 
	pos = sram$pos
	repeats = sram$repeatname
	rnames = sram$rname

	cb = which(diff(pos,1) > gap.cutoff)
	cb = c(cb, length(pos))

	print(paste(length(cb), "raw breakpoints"))
	cl = data.frame(do.call(rbind, lapply(1:length(cb), function(i) {

		if (verbose & i%%100 == 0) print(paste("procissing", i))
 
		if (i==1) si=1 else si = cb[i-1] + 1
		ei = cb[i]
		rv = unlist(lapply(repeats[si:ei], function(r) {
			if (length(grep(",", r))) return(strsplit(r, ",")[[1]][[1]]) else return(r)
		}))
		rep.repeat = paste(unique(rv), collapse=",")

		if (merge.family) {
			rannot = load.rannot()
	  		rfamily = rannot$family[match(toupper(rv), rannot$name)]; rfamily[is.na(rfamily)] = "NA"
	  		rclass = rannot$class[match(toupper(rv), rannot$name)]; rclass[is.na(rclass)] = "NA"
		} else {
			rfamily = rv
			rclass = rep("-", length(rfamily))
		}

	  	cl.row = data.frame(s=pos[si], e=pos[ei], ram=ei-si+1, rep.repeat=rep.repeat, family=paste(rfamily, collapse=","), class=paste(rclass, collapse=","), repeat.name=paste(rv, collapse=","), pos=paste(pos[si:ei]*strand, collapse=","), rname=paste(rnames[si:ei], collapse=","), stringsAsFactors=F)

		separated.cl.row = separate.cluster2(as.list(cl.row))

	  	return(separated.cl.row)
	})))
	cl$s = as.integer(as.character(cl$s))
	cl$e = as.integer(as.character(cl$e))
	cl$ram = as.integer(as.character(cl$ram))
	cl$rep.repeat = as.character(cl$rep.repeat)
	cl$family = as.character(cl$family)
	cl$class = as.character(cl$class)
	cl$rname = as.character(cl$rname)

	return(cl)
}

# change the PolyA to the majority among L1, Alu, SVA when they are mixed
get.cluster3 <- function(sram, strand, gap.cutoff, merge.family=F, verbose=F)
{ 
	pos = sram$pos
	repeats = sram$repeatname
	rnames = sram$rname

	cb = which(diff(pos,1) > gap.cutoff)
	cb = c(cb, length(pos))

	print(paste(length(cb), "raw breakpoints"))
	cl = data.frame(do.call(rbind, lapply(1:length(cb), function(i) {

		if (verbose & i%%100 == 0) print(paste("procissing", i))
 
		if (i==1) si=1 else si = cb[i-1] + 1
		ei = cb[i]
		rv = unlist(lapply(repeats[si:ei], function(r) {
			if (length(grep(",", r))) return(strsplit(r, ",")[[1]][[1]]) else return(r)
		}))
		rep.repeat = paste(unique(rv), collapse=",")

		if (merge.family) {
			rannot = load.rannot()
	  		rfamily = rannot$family[match(toupper(rv), rannot$name)]; rfamily[is.na(rfamily)] = "NA"
	  		rclass = rannot$class[match(toupper(rv), rannot$name)]; rclass[is.na(rclass)] = "NA"

   			# convert the family of PolyA to the majority among L1, Alu, or SVA when they're mixed
            if (grepl("PolyA", rv) && length(unique(rv))>1) {
                # identify the majority family among L1, Alu and SVA
				x = rfamily[rfamily == "L1" | rfamily == "Alu" | rfamily == "SVA"]
                newfamily = names(table(x))[table(x) == max(table(x))][1]    
                rv[rv=="PolyA"] = newfamily
            }
	  		rfamily = rannot$family[match(toupper(rv), rannot$name)]; rfamily[is.na(rfamily)] = "NA"
	  		rclass = rannot$class[match(toupper(rv), rannot$name)]; rclass[is.na(rclass)] = "NA"
		} else {
			rfamily = rv
			rclass = rep("-", length(rfamily))
		}

	  	cl.row = data.frame(s=pos[si], e=pos[ei], ram=ei-si+1, rep.repeat=rep.repeat, family=paste(rfamily, collapse=","), class=paste(rclass, collapse=","), repeat.name=paste(rv, collapse=","), pos=paste(pos[si:ei]*strand, collapse=","), rname=paste(rnames[si:ei], collapse=","), stringsAsFactors=F)

		separated.cl.row = separate.cluster2(as.list(cl.row))

	  	return(separated.cl.row)
	})))
	cl$s = as.integer(as.character(cl$s))
	cl$e = as.integer(as.character(cl$e))
	cl$ram = as.integer(as.character(cl$ram))
	cl$rep.repeat = as.character(cl$rep.repeat)
	cl$family = as.character(cl$family)
	cl$class = as.character(cl$class)
	cl$rname = as.character(cl$rname)

	return(cl)
}
separate.cluster2 <- function(cl.row)
{
	uq.rfamily = unique(strsplit(cl.row$family, ",")[[1]])

	# genearte rows for each cluster of one family
	separated.cl.row = do.call(rbind, lapply(1:length(uq.rfamily), function(i) {
		idx = which(strsplit(cl.row$family, ",")[[1]] == uq.rfamily[i])
	  	pos = strsplit(cl.row$pos, ",")[[1]][idx]
	  	repeat.name = strsplit(cl.row$repeat.name, ",")[[1]][idx]
		rname = strsplit(cl.row$rname, ",")[[1]][idx]
		rfamily = uq.rfamily[i] 
		rclass = unique(strsplit(cl.row$class, ",")[[1]][idx])
		s = abs(as.numeric(head(pos, 1))); e = abs(as.numeric(tail(pos, 1)))
	  	ram = length(pos)
	  	rep.repeat=unique(repeat.name)
	  	return(c(s=s, e=e, ram=ram, rep.repeat=paste(rep.repeat, collapse=","),  family=paste(rfamily, collapse=","), class=paste(rclass, collapse=","), repeat.name=paste(repeat.name, collapse=","), pos=paste(pos, collapse=","), rname=paste(rname, collapse=",")))
	}))
	return(separated.cl.row)
}

pair.cluster2 <- function(p.cl, m.cl, gap.cutoff, stringent.pair=F, verbose=F)
{
	ov.raw = as.matrix(findOverlaps(IRanges(p.cl$s, p.cl$e + gap.cutoff), IRanges(m.cl$s, m.cl$e)))
	colnames(p.cl) = paste(colnames(p.cl), "1", sep="")
	colnames(m.cl) = paste(colnames(m.cl), "2", sep="")
	pcl.raw = cbind(p.cl[ov.raw[,1],], m.cl[ov.raw[,2],])

	# exclude inappropriate pairs such as
    # if the negative strand cluster appear on the left side of the positive stand cluster
    # if family is annotated as "Unknown" 
	if (stringent.pair) {
		idx = which(pcl.raw$s1 <= pcl.raw$s2 & pcl.raw$e1 <= pcl.raw$e2 & pcl.raw$family1 == pcl.raw$family2 & 
			pcl.raw$family1 != "Unknown" & pcl.raw$family2 != "Unknown")
	} else {
		idx = which(pcl.raw$s1 < pcl.raw$e2 & pcl.raw$family1 == pcl.raw$family2 & 
			pcl.raw$family1 != "Unknown" & pcl.raw$family2 != "Unknown")
	}
	ov = ov.raw[idx, ,drop=F]
	pcl = pcl.raw[idx, ,drop=F]

	return(list(ov=ov, pcl=pcl))
}

# consider reads mapped to polyA to be either L1, Alu, or SVA family
pair.cluster3 <- function(p.cl, m.cl, gap.cutoff, stringent.pair=F, verbose=F)
{
    ov.raw = as.matrix(findOverlaps(IRanges(p.cl$s, p.cl$e + gap.cutoff), IRanges(m.cl$s, m.cl$e)))
    colnames(p.cl) = paste(colnames(p.cl), "1", sep="")
    colnames(m.cl) = paste(colnames(m.cl), "2", sep="")
    pcl.raw = cbind(p.cl[ov.raw[,1],], m.cl[ov.raw[,2],])

    # exclude inappropriate pairs such as
    # if the negative strand cluster appear on the left side of the positive stand cluster
    # if family is annotated as "Unknown" 
    if (stringent.pair) {
        idx = which(pcl.raw$s1 <= pcl.raw$s2 & pcl.raw$e1 <= pcl.raw$e2 &
            (pcl.raw$family1 == "PolyA" | pcl.raw$family2 == "PolyA" |
            pcl.raw$family1 == pcl.raw$family2) &
            pcl.raw$family1 != "Unknown" & pcl.raw$family2 != "Unknown")
    } else {
        idx = which(pcl.raw$s1 < pcl.raw$e2 &
            (pcl.raw$family1 == "PolyA" | pcl.raw$family2 == "PolyA" |
            pcl.raw$family1 == pcl.raw$family2) &
            pcl.raw$family1 != "Unknown" & pcl.raw$family2 != "Unknown")
    }
    ov = ov.raw[idx, ,drop=F]
    pcl = pcl.raw[idx, ,drop=F]

    return(list(ov=ov, pcl=pcl))
}

get.cluster <- function(pos, strand, repeats, rnames, gap.cutoff, family.input=F, verbose=F)
{ 
	cb = which(diff(pos,1) > gap.cutoff)
	cb = c(cb, length(pos))

	print(paste(length(cb), "raw breakpoints"))
	cl = data.frame(do.call(rbind, lapply(1:length(cb), function(i) {

		if (verbose & i%%100 == 0) print(paste("procissing", i))
 
		if (i==1) si=1 else si = cb[i-1] + 1
		ei = cb[i]
		rv = unlist(lapply(repeats[si:ei], function(r) {
			if (length(grep(",", r))) return(strsplit(r, ",")[[1]][[1]]) else return(r)
		}))
		rep.repeat = paste(unique(rv), collapse=",")

		if (!family.input) {
			rannot = load.rannot()
	  		rfamily = rannot$family[match(rv, rannot$name)]; rfamily[is.na(rfamily)] = "NA"
	  		rclass = rannot$class[match(rv, rannot$name)]; rclass[is.na(rclass)] = "NA"
		} else {
			#debug end
			#rfamily = strsplit(rep.repeat, ",")[[1]]
			#rclass = rfamily 
			rfamily = rv
			rclass = rep("-", length(rfamily))
			#debut end
		}

	  	cl.row = c(s=pos[si], e=pos[ei], ram=ei-si+1, rep.repeat=rep.repeat, family=paste(rfamily, collapse=","), 
			class=paste(rclass, collapse=","), repeat.name=paste(rv, collapse=","), pos=paste(pos[si:ei]*strand, 
			collapse=","), rname=paste(rnames[si:ei], collapse=","))

		# separate a cluster if it's mixed with the repeats of the different repeat family
		if (length(unique(rfamily))>1) {
			cl.row = separate.cluster(cl.row, rfamily)
		} else {
			cl.row = c(s=pos[si], e=pos[ei], ram=ei-si+1, rep.repeat=rep.repeat, family=unique(rfamily), 
			class=unique(rclass), repeat.name=paste(rv, collapse=","), pos=paste(pos[si:ei]*strand, collapse=","), 
			rname=paste(rnames[si:ei], collapse=","))
	  	}
		print(as.list(cl.row)$rep.repeat)
		print(as.list(cl.row)$family)
	  	return(cl.row)
	})))
	cl$s = as.integer(as.character(cl$s))
	cl$e = as.integer(as.character(cl$e))
	cl$ram = as.integer(as.character(cl$ram))
	cl$rep.repeat = as.character(cl$rep.repeat)
	cl$family = as.character(cl$family)
	cl$class = as.character(cl$class)
	cl$rname = as.character(cl$rname)

	return(cl)
}

pair.cluster <- function(p.cl, m.cl, gap.cutoff)
{
	ov.raw = as.matrix(findOverlaps(IRanges(p.cl$s, p.cl$e + gap.cutoff), IRanges(m.cl$s, m.cl$e)))
	colnames(p.cl) = paste(colnames(p.cl), "1", sep="")
	colnames(m.cl) = paste(colnames(m.cl), "2", sep="")
	pcl.raw = cbind(p.cl[ov.raw[,1],], m.cl[ov.raw[,2],])

	# exclude pairs if the negative strand clusters appear on the left side of the positive stand clusters: s1 < e2
	idx = which(pcl.raw$s1 < pcl.raw$e2 & pcl.raw$family1 == pcl.raw$family2 & pcl.raw$family1 != "Unknown" 
		& pcl.raw$family2 != "Unknown")
	ov = ov.raw[idx, ,drop=F]
	pcl = pcl.raw[idx, ,drop=F]
	# check the repeat family are the same for the positive and negative clusters; but dont pair if the family is "unknown"
	#ovidx = unlist(lapply(1:dim(pcl.raw)[1], function(i) { 
	# p.rf = pcl.raw$family1[i]
	# m.rf = pcl.raw$family2[i]
	# if (p.rf != "Unknown" & m.rf != "Unknown" & p.rf == m.rf) return(i) else return(NULL)
	#}))
	#pcl = pcl.raw[ovidx,]

	return(list(ov=ov, pcl=pcl))
}

separate.cluster <- function(cl.row, rf)
{
	cl.row = as.list(cl.row)
	unique.rf = unique(rf)

	separated.cl.row = do.call(rbind, lapply(1:length(unique.rf), function(i) {
	  idx = which(rf == unique.rf[i])
	  pos = strsplit(cl.row$pos, ",")[[1]][idx]
	  repeat.name = strsplit(cl.row$repeat.name, ",")[[1]][idx]
	  rname = strsplit(cl.row$rname, ",")[[1]][idx]
	  rep.repeat=unique(repeat.name)
	  rfamily = strsplit(cl.row$family, ",")[[1]][idx]
	  rclass = strsplit(cl.row$class, ",")[[1]][idx]
	  s = abs(as.numeric(head(pos, 1))); e = abs(as.numeric(tail(pos, 1)))
	  ram = length(pos)
	  return(c(s=s, e=e, ram=ram, rep.repeat=paste(rep.repeat, collapse=","),  family=paste(unique(rfamily), collapse=","), class=paste(unique(rclass), collapse=","), repeat.name=paste(repeat.name, collapse=","), pos=paste(pos, collapse=","), rname=paste(rname, collapse=",")))
	}))
	return(separated.cl.row)
}

load.rannot <- function(rannot.file=NULL)
{
	if (is.null(rannot.file)) rannot.file = sprintf("%s/lib/repeats.txt", tea.base)

	rannot = read.table(rannot.file, sep="\t", header=F, col.names=c("name","class"), as.is=T)
	rannot$class.org = rannot$class
	rannot$class = unlist(lapply(rannot$class.org, function(x) strsplit(x, "/")[[1]][[1]]))
	rannot$family = unlist(lapply(rannot$class.org, function(x) {
	    if (length(strsplit(x, "/")[[1]]) > 1) {
	      return(strsplit(x, "/")[[1]][[2]])
	    } else  {
	      return("NA")
	    }
	    }))
	return(rannot)
}

get.score <-function(can, max.ram=20, max.cr=10)
{
    s.r1 = mapply(min, as.numeric(can$pram), max.ram)/max.ram
    s.r2 = mapply(min, as.numeric(can$nram), max.ram)/max.ram

	if (("pacr" %in% colnames(can))) { 
    	s.c1 = mapply(min, as.numeric(can$pacr), max.cr)/max.cr
    	s.c2 = mapply(min, as.numeric(can$nacr), max.cr)/max.cr
	} else { # for column name compatibility
		s.c1 = mapply(min, as.numeric(can$acr1), max.cr)/max.cr
		s.c2 = mapply(min, as.numeric(can$acr2), max.cr)/max.cr
	}

    score = s.r1 + s.r2 + s.c1 + s.c2

    return(score)
}
 
# return two data frames for cluster (with pos) and cluster detail (with pos and tname)
cluster.ram <- function(desc, ram, is, ril, rname, rfamily, min.single.ram=3, verbose=F)
{
	chrl = names(ram); names(chrl) = chrl

	cl = data.frame(do.call(rbind, lapply(chrl, function(c) {
	#cl = data.frame(do.call(rbind, lapply(chrl[1:2], function(c) {
	  if (verbose) print(paste("processing", c))
	  if (dim(ram[[c]]$p)[1] == 0 | dim(ram[[c]]$m)[1] == 0) return (NULL)
	  p = data.frame(rname=ram[[c]]$p$rname, strand=1, pos=ram[[c]]$p$pos)
	  m = data.frame(rname=ram[[c]]$m$rname, strand=-1, pos=ram[[c]]$m$pos)
	  pm = rbind(p,m)
	  pm = pm[order(pm$pos),]
	  tram = dim(pm)[1]

	  # identifying cluster boundaries
	  cb = c(which((pm$pos[-1] - pm$pos[-tram]) > is$gap.size), tram)
	  if (verbose) write.msg(paste(length(cb), "cluster boundaries"))

	  insf=NULL
	  ins = ril[[rname]][[c]]
	  if (!is.null(rfamily)) {
	    insf = do.call(rbind, lapply(rfamily, function(r) {
	      return(ril[[r]][[c]])
	    }))
	  }

	  df = data.frame(do.call(rbind, lapply(1:length(cb), function(i) {
	    #if (verbose) print(paste("processing", i))
	    si = ifelse(i==1, 1, cb[i-1]+1)
	    ei = cb[i]
	    ram = ei - si + 1
	    pram = length(which(pm$strand[si:ei] >0))
	    nram = length(which(pm$strand[si:ei] <0))
	    s = pm$pos[si]
	    e = pm$pos[ei]
	    size = e-s+1
	    return(c(chr=c, s=s, e=e, size=size, ram=ram, pram=pram, nram=nram))
	  })))

	  if (dim(df)[1] != 0) {
	    df$s = as.integer(as.character(df$s))
	    df$e = as.integer(as.character(df$e))
	    df$size = as.integer(as.character(df$size))
	    df$ram = as.integer(as.character(df$ram))
	    df$pram = as.integer(as.character(df$pram))
	    df$nram = as.integer(as.character(df$nram))

	    # adjusting cluster boundaries based on the peak patterns
	    df$s[df$pram == 0] = df$s[df$pram == 0] - (is$mu + 3 * is$id)
	    df$e[df$nram == 0] = df$e[df$nram == 0] + desc$rl + (is$mu + 3 * is$sd)
	    df$size = df$e - df$s

	    # known instance annotaton
	    df$oi = ifelse(countOverlaps(IRanges(df$s, df$e), IRanges(ins$s-is$ins.margin, ins$e+is$ins.margin), type="any")==0,1,0)
	    df$oif = ifelse(countOverlaps(IRanges(df$s, df$e), IRanges(insf$s-is$ins.margin, insf$e+is$ins.margin), type="any")==0,1,0)

	    # define regions for clipped read search 
	    # cs, ce: clipped search start,end 
	    # ls, rs: left/right cluster start
	    # le, re: left/right cluster end 
	    df$cs = -1; df$ce = -1;  df$ls = -1; df$le = -1; df$rs = -1; df$re = -1
	    df$pos = ""; df$tname = ""

	    si = head(c(1, cb)+1, -1)
	    si[1] = 1
	    ei = cb
	    idx = which(df$oi == 1 & (df$pram >= min.single.ram | df$nram >=min.single.ram) )

	    if (verbose) write.msg(paste(rname, c, length(df$s), "total", length(idx), "oi clusters having at least", min.single.ram, "rams"))

	    if (length(idx) >= 1) {
	      for (i in 1:length(idx)) {
	        #print(paste("processing", i))
	        spm = pm[si[idx[i]]:ei[idx[i]],]
	        df$tname[idx[i]] = paste(as.character(spm$rname), collapse=" ")
	        df$pos[idx[i]] = paste(spm$pos * spm$strand, collapse=" ")

	        if (df$pram[idx[i]] >0 & df$nram[idx[i]] >0) {
	            bp = get.bp("both", spm$pos*spm$strand, desc$rl, wrl=100, max.tsd=50, margin=50)
	        } else if (df$nram[idx[i]] == 0) {
	           bp = get.bp("left", spm$pos*spm$strand, desc$rl, wrl=100, max.tsd=50, margin=50)
	        } else {
	            bp = get.bp("right", spm$pos*spm$strand, desc$rl, wrl=100, max.tsd=50, margin=50)
	        }
	        df[idx[i], c("cs", "ce", "ls", "le", "rs", "re")] = bp
	      }
	    }
	    rm(p, m, pm, ins)
	    return(df)
	  } else {
	    rm(p, m, pm, ins)
	    return(NULL)
	  }
	})))

	if (dim(cl)[1] != 0) {
	    return(cl)
	} else {
	  return (NULL)
	}
}

count.clipped2 <- function(can, bamf, verbose=F, jittering=0, cbam.chr=T, rl=100)
{
	if (!cbam.chr)	can$chr = gsub("chr", "", can$chr)
	
	start = can$pram_end
	start[can$pram_end == 0] = Inf
	start = pmin(can$s + rl, start - rl)
	end = pmax(can$nram_start + rl, can$e - rl)

	what <- c("qname", "pos", "strand", "cigar", "seq", "qual")

	clipped = lapply(1:length(can$s), function(i) {
	#clipped = lapply(1:20, function(i) {
		if (verbose) print(paste("processing", i))
		#which<- GRanges(seqnames = can$chr[i], ranges = IRanges(can$s[i], can$e[i])-rl)
		which <- GRanges(seqnames = can$chr[i], ranges = IRanges(start[i], end[i]))
		param <- ScanBamParam(which=which, what=what)
		map <- scanBam(bamf, param=param)

		ccnt = list(acr=0, acr1=0, acr2=0, cr=0, cr1=0, cr2=0, pbp=-1, nbp=-1, tsd=-9999, arr=0)
		cdf= NULL
		if (length(map[[1]]$seq)>0 & length(map[[1]]$seq)<1000) {
			lst <- lapply(names(map[[1]]), function(elt) {
		    do.call("c", unname(lapply(map, "[[", elt)))
			})
			names(lst) <- names(map[[1]])
			map.df <- do.call("DataFrame", lst)

			clipped.df1 = get.clipped.pos(map.df, qual.trim=T)
			clipped.df2 = get.clipped.pos(map.df, qual.trim=F)

			if (nrow(clipped.df1) > 0)
			ccnt1 = get.clipped.cnt2(as.integer(as.character(clipped.df1$cpos)), jittering)

			if (nrow(clipped.df2) > 0)
			ccnt2 = get.clipped.cnt2(as.integer(as.character(clipped.df2$cpos)), jittering)

			if (ccnt1$acr >= ccnt2$acr)
				realigned.clipped = realign.clipped(ccnt1, clipped.df1, clipped.df2, jittering)
		  	else
		    	realigned.clipped = realign.clipped(ccnt2, clipped.df1, clipped.df2, jittering)

		  	ccnt = realigned.clipped$ccnt
			cdf = merge(can[i, match(c("chr", "s", "e"), colnames(can))],realigned.clipped$cdf)
		}
		if (verbose) print(paste("processing", i, ccnt$cr, ccnt$cr1, ccnt$cr2, ccnt$acr, ccnt$acr1, ccnt$acr2, ccnt$arr))
		return(list(ccnt=ccnt, cdf=cdf))
	})
	return(clipped)
}

realign.clipped <- function(ccnt, df1, df2, jittering)
{
	df1$aligned = 0; df2$aligned = 0
	cdf = data.frame(do.call(rbind, lapply(1:dim(df1)[1], function(i) {
	  if (df1$cpos[i] > 0) {
	    if (abs(df1$cpos[i] - ccnt$pbp) <= jittering) {
	      df1$aligned[i] = 1
	      return(df1[i,])
	    } else if (abs(df2$cpos[i] - ccnt$pbp) <= jittering) {
	      df2$aligned[i] = 1
	      return(df2[i,])
	    } else {
	      return(df1[i,])
	    }
	  } else {
	    if (abs(df1$cpos[i] + ccnt$nbp) <= jittering) {
	      df1$aligned[i] = 1
	      return(df1[i,])
	    } else if (abs(df2$cpos[i] + ccnt$nbp) <= jittering) {
	      df2$aligned[i] = 1
	      return(df2[i,])
	    } else {
	      return(df1[i,])
	    }
	  }
	})))
	ccnt$acr1 = length(which(cdf$aligned == 1 & cdf$cpos >0))
	ccnt$acr2 = length(which(cdf$aligned == 1 & cdf$cpos <0))
	ccnt$acr = ccnt$acr1 + ccnt$acr2
	return(list(ccnt=ccnt, cdf=cdf))
}

# return maxcnt and clipped position
get.clipped.cnt3 <- function(pos, jittering)
{
	cnt=length(pos); acnt=0; cp=-1
	if (cnt==0) return(list(cnt=cnt, acnt=acnt, cp=cp))

	x = pos
	for (i in 1:jittering) x = c(x, pos-i, pos+i)

	acnt = max(table(x))

	org.maxidx = which(table(pos) == max(table(pos)))
	maxidx = which(table(x) == max(table(x)))

	#if (names(org.maxidx) %in% names(maxidx)) {
	  if (length(org.maxidx) ==1) cp = as.integer(names(org.maxidx)) else
	    cp = as.integer(names(org.maxidx)[floor(length(org.maxidx)/2)])
	#} else {
	# cp = as.integer(names(maxidx)[floor(length(maxidx)/2)]) 
	#}
	return(list(cnt=cnt, acnt=acnt, cp=cp))
}

get.clipped.cnt2 <- function(cpos, jittering=0)
{
	p = get.clipped.cnt3(cpos[cpos>0], jittering)
	m = get.clipped.cnt3(-cpos[cpos<0], jittering)

	cr1 = p$cnt; cr2 = m$cnt; cr = cr1+cr2
	acr1 = p$acnt; acr2 = m$acnt; acr = acr1 + acr2
	pbp = p$cp; nbp = m$cp
	arr = round(acr / cr, 2)
	tsd = -9999
	if (nbp != -1 & pbp != -1)  tsd = nbp - pbp  - 1

	return(list(acr = acr, acr1 = acr1, acr2 = acr2, cr=cr, cr1=cr1, cr2=cr2, pbp=pbp, nbp=nbp, tsd=tsd, arr=arr))
}

get.clipped.pos <- function(map.df, qual.trim=F)
{

	df = data.frame(do.call(rbind, lapply(1:dim(map.df)[1], function(i) {
	  d = map.df[i,]
	  rclip.len=0; lclip.len=0

	  if (length(grep(".*[MIDNHPS](\\d+)S$", d$cigar))!=0)
	    rclip.len = as.integer(sub(".*[MIDNHPS](\\d+)S$", "\\1", d$cigar, perl=T))
	  if (length(grep("^(\\d+)S.*", d$cigar))!=0)
	    lclip.len = as.integer(sub("^(\\d+)S.*", "\\1", d$cigar, perl=T))

	  #print(paste(i, d$cigar, lclip.len, rclip.len))
	  if (d$strand == 1 & lclip.len >0) {
	    cpos = d$pos - 1
	    ref.seq = substr(d$seq[[1]], lclip.len+1, length(d$seq[[1]])-rclip.len)
	    clipped.seq = substr(d$seq[[1]], 1, lclip.len)
	    clipped.qual = substr(d$qual[[1]], 1, lclip.len)

	    if (qual.trim & length(grep("[^#]*(#+)$", clipped.qual)) != 0) {
	      # count # on the right of qual
	      badqual.cnt = nchar(sub("[^#]*(#+)$", "\\1", clipped.qual, perl=T))
	      cpos = cpos - badqual.cnt
	      clipped.seq = substr(clipped.seq, 1, nchar(clipped.seq) - badqual.cnt)
	      clipped.qual = substr(clipped.qual, 1, nchar(clipped.qual) - badqual.cnt)
	    }

	    return(c(cpos=cpos, cigar=d$cigar, tname=d$qname, ref.seq=toString(ref.seq), clipped.seq=toString(clipped.seq), clipped.qual = toString(clipped.qual)))
	  } else if (d$strand ==2 & rclip.len >0) {
	    # count deletion
	    del.cnt = 0; new.cigar = d$cigar
	    while (length(grep(".*[MINHPS](\\d+)D.*", new.cigar)) != 0) {
	      del.cnt = del.cnt + as.integer(sub(".*[MINHPS](\\d+)D.*", "\\1", new.cigar))
	      new.cigar = sub("(.*)\\d+D.*", "\\1", new.cigar)
	    }
	    ins.cnt = 0; new.cigar = d$cigar
	    while (length(grep(".*[MDNHPS](\\d+)I.*", new.cigar)) != 0) {
	      ins.cnt = ins.cnt + as.integer(sub(".*[MDNHPS](\\d+)I.*", "\\1", new.cigar))
	      new.cigar = sub("(.*)\\d+I.*", "\\1", new.cigar)
	    }
	    cpos = -(as.integer(d$pos) + length(d$seq[[1]]) - lclip.len - rclip.len + del.cnt - ins.cnt)
	    ref.seq = substr(d$seq[[1]], lclip.len+1, length(d$seq[[1]]) - rclip.len)
	    clipped.seq = substr(d$seq[[1]], length(d$seq[[1]]) - rclip.len + 1, length(d$seq[[1]]))
	    clipped.qual = substr(d$qual[[1]], length(d$qual[[1]]) - rclip.len + 1, length(d$qual[[1]]))

	    if (qual.trim & length(grep("^(#+)[^#]*", clipped.qual)) != 0) {
	      # count # on the left of the qual
	      badqual.cnt = nchar(sub("^(#+)[^#]*", "\\1", clipped.qual, perl=T))
	      cpos = cpos - badqual.cnt
	      clipped.seq = substr(clipped.seq, badqual.cnt+1, length(clipped.seq))
	      clipped.qual = substr(clipped.qual, badqual.cnt+1, length(clipped.qual))
	    }

	    return(c(cpos=cpos, cigar=d$cigar, tname=d$qname, ref.seq=toString(ref.seq), clipped.seq=toString(clipped.seq), clipped.qual=toString(clipped.qual)))
	  } else {
	    return(NULL)
	  }
	})))
	df$cpos = as.integer(as.character(df$cpos))
	return(df)
}


load.rl <- function(rl.file, rg="all")
{
	if (is.null(rl.file)) return(NULL)
	rl = read.table(rl.file)
	return(rl[rl$V1=="all", 2])
}

load.isize <- function(isize.file=NULL, rl, mu, sd)
{
	if (!is.null(isize.file)) {
	  is = read.table(isize.file)
	  mu = floor(is[is$V1=="all",2])
	  sd = is[is$V1=="all",3]
	}
	fr = floor(mu + 3 * sd)
	intra.gap = floor((fr-rl) * 2/3)
	inter.gap = intra.gap * 2
	gap.size = floor(mu*2/3 + rl)
	ins.margin = floor(1.5 * fr)
	return(list(fr = fr, mu=mu, sd=sd, gap.size=gap.size, intra.gap=intra.gap, inter.gap=inter.gap, ins.margin=ins.margin))
}

load.ref.info <- function(ref="hg18", out.chrl=T, out.ril=F, out.gene=F, out.gap=F, merged.rmasker=T, verbose=F)
{
	chrl = ril = gap = genes = NULL
	if (verbose) write.msg(paste("loading information for the reference", ref, "..."))

	# chromosome list 
	if (out.chrl) {
		if (ref == "hg18" | ref == "hg19") chrl = paste("chr", c(1:22, "X", "Y"), sep="") else
		if (ref == "hg18" | ref == "hg19") chrl = paste("chr", c(1:22, "X", "Y"), sep="") else
		if (ref == "ponAbe2") chrl = paste("chr", c(1, "2a", "2b", 3:22, "X"), sep="") else
		if (ref == "panTro3") chrl = paste("chr", c(1, "2A", "2B", 3:22, "X", "Y"), sep="") else
		if (ref == "rheMac2") chrl = paste("chr", c(1:20, "X"), sep="")
		names(chrl) = chrl
		if (verbose) write.msg(paste(length(chrl), "chrs:", paste(chrl, collapse=",")))
	}

	# reference te instances annotated by repeatmasker
	if (out.ril) {
		if (merged.rmasker) {
			rmasker.rfile = sprintf("%s/lib/rmasker/rmasker.%s.merged.RData", tea.base, ref)
		} else {
			rmasker.rfile = sprintf("%s/lib/rmasker/rmasker.%s.RData", tea.base, ref)
		}
		if (file.exists(rmasker.rfile)) {
			print(rmasker.rfile)
			write.msg(sprintf("loading rmasker te annotation %s", rmasker.rfile))
			ril = readRDS(rmasker.rfile) 
			write.msg("done.")
		} else {
			stop(sprintf("no TE instance annoation file: %s", rmasker.rfile))
		}
		write.msg(sprintf("known repeatmasker te instances were loaded for %d repeat types", length(ril)))
	}
	
	# reference gap annotation
	if (out.gap) {
		gap.fname = sprintf("%s/lib/gap/gap.%s.gz", tea.base, ref)
		write.msg(sprintf("reading %s .......", gap.fname))
		gap = read.gap.rfile(gap.fname, chrl)
		write.msg(sprintf("done %d gapped regions have been loaded.\n", sum(sapply(gap, nrow))))
	}

	# gene annotation
	if (out.gene) {
		gene.rfile = sprintf("%s/lib/gene/%s.genes.RData", tea.base, ref)
		if (file.exists(gene.rfile)) {
			write.msg(sprintf("loading the gene annotation %s.......", gene.rfile))
			genes = readRDS(gene.rfile)
			write.msg("done.")
		} else {
			stop(sprintf("no gene annoation file: %s", gene.rfile))
		}
	}
	ref = list(chrl = chrl, ril=ril, gap=gap, genes = genes)

	return(ref)
}

# input repeatmasker files were generated using "Tables" in the UCSC genome browser web site
make.rmasker.rfile <- function(rmasker.rfile, chrl=NULL)
{
	rmasker.fname = sub(".RData", ".gz", rmasker.rfile)
	coordinate.fname = sub(".RData", ".coordinates.gz", rmasker.rfile)
	if (!file.exists(coordinate.fname)) {
		# generate te instance coordinates using rmasker_coordinates (ra.pl)
		cmd = sprintf("%s/scripts/ra.pl rmasker_coordinates %s 2", tea.base, rmasker.fname)
		write.msg("generating te instance annotation.......")
		write.msg(sprintf("%s", cmd))
		system(cmd)
		write.msg("done\n")
	}

	 rt <- read.delim(gzfile(coordinate.fname), sep="\t", header=F, as.is=T)
	 ril <- lapply(1:dim(rt)[1],function(i) rt[i,2]); names(ril) <- rt[,1]
	 ril <- lapply(ril, function(s) {
	   d <- do.call(rbind,lapply(strsplit(s," ")[[1]],function(x) strsplit(x,":|-")[[1]]))
			if (!is.null(chrl)) d = d[d[,1] %in% chrl, ,drop=F ]
			#print(paste(strsplit(s," ")[[1]][[1]]))
			if (nrow(d) == 0) return(NULL)
	 		return(tapply(1:dim(d)[1],as.factor(d[,1]), function(ii)
	             data.frame(s=as.integer(d[ii,2]),e=as.integer(d[ii,3]))))
	 })

	write.msg(sprintf("creating rmasker.rfile %s........", rmasker.rfile))
	saveRDS(ril, file=rmasker.rfile)
	write.msg(sprintf("done: known repeatmasker te instances were loaded for %d repeat types\n", length(ril)))


	cmd = sprintf("rm %s", coordinate.fname)
	system(cmd)
	write.msg(sprintf("cleaning %s done", coordinate.fname))

	return(ril)
}

# input gap files were generated using "Tables" in the UCSC genome browser web site 
read.gap.rfile <- function(gap.fname, chrl=NULL)
{
	gap = read.delim(gzfile(gap.fname), sep="\t", as.is=T)
	if (grepl("rheMac2", gap.fname)) gap = gap[, 1:3] else gap = gap[, 2:4]
	colnames(gap) = c("chr", "s", "e")
	gap$s = gap$s + 1 # note that the start coordiates are 0-based
	
	if (!is.null(chrl)) gap = gap[gap$chr %in% chrl, ]

 	gap = tapply(1:dim(gap)[1],as.factor(gap[,1]), function(ii)
	             data.frame(s=as.integer(gap[ii,2]),e=as.integer(gap[ii,3])))

	return(gap)
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

write.msg <- function(msg)
{
	write(msg, stdout())
}


# calculate signal to noise ratio 
# 2*sqrt(#pram_downstream_of_nbp * #nram_upstream_of_pbp) - (#pram_upstreadm_of_nbp + #nram_downstream_of_pbp)
# assume cl.file is only for one sample
s2n <- function(dir, sample, rasym="ra", cl.file=NULL, ref="hg19")
{
    ref.annot = load.ref.info(ref)
    chrl = ref.annot$chrl

    # load cluster and raw_cluster files 
    cldir = sprintf("%s/%s/cluster_%sm", dir, sample, rasym)
    cl.df = read.delim(paste(cldir, "/", cl.file, sep=""), sep="\t", as.is=T)
    cl = split(cl.df, cl.df$chr)
    out.file = sprintf("%s/%s.s2n", cldir, cl.file)

    # for each chr, calculate s2n
    x = do.call(rbind, lapply(names(cl), function(chr) {
        raw.file = sprintf("%s/%s/cluster_%sm/%s.%s.cluster.raw", dir, sample, rasym, sample, chr)
        raw = read.delim(raw.file, header=T, as.is=T)
        m = merge(cl[[chr]], raw, by=c("s", "e"))
        m = m[, c("pbp", "nbp", "pos1", "pos2")]

        cl[[chr]]$s2n = unlist(lapply(1:nrow(m), function(i) calc.s2n(m[i,])))
        return(cl[[chr]])
    }))
    write.table(x, out.file, sep="\t", quote=F, row.names=F)
    return(x)
}

calc.s2n <- function(x) {
	pram = as.numeric(strsplit(x$pos1, ",")[[1]])
	nram = -as.numeric(strsplit(x$pos2, ",")[[1]])

	pram_down_nbp = length(which(pram <= x$nbp))
	pram_up_nbp = length(which(pram > x$nbp))

	nram_up_pbp = length(which(nram >= x$pbp))
	nram_down_pbp = length(which(nram < x$pbp))

	s2n = 2*sqrt(pram_down_nbp * nram_up_pbp) - (pram_up_nbp + nram_down_pbp)
	return(s2n)
}
