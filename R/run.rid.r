tea.base = Sys.getenv("tea_base")
if (tea.base == "")  stop("tea.base environment variable needs to be set") 
source(paste(tea.base, "/R/rid.r", sep=""))
args <- commandArgs(trailingOnly=T)

if (args[1] == "rid") {

	chr=NULL; no.clipped=F; ref="hg18"; rasym="ra"; min.ram=3; oneside.ram=F; exo=F; cbam.chr=F
	jittering=2; merged.rmasker=T; merge.family=T; annot.oi=T; annot.gene=T; rmasker.filter.margin=500; 
	verbose=F

	for(i in 2:length(args)) eval(parse(text=args[[i]]))

	write.msg(paste("dir:", dir))
	write.msg(paste("sample:", sample))
	write.msg(paste("chr:", chr))
	write.msg(paste("no.clipped:", no.clipped))
	write.msg(paste("ref:", ref))
	write.msg(paste("rasym:", rasym))
	write.msg(paste("min.ram:", min.ram))
	write.msg(paste("oneside.ram:", oneside.ram))
	write.msg(paste("exo:", exo))
	write.msg(paste("cbam.chr:", F))
	write.msg(paste("jittering:", jittering))
	write.msg(paste("merged.rmasker:", merged.rmasker))
	write.msg(paste("merge.family:", merge.family))
	write.msg(paste("annot.oi:", annot.oi))
	write.msg(paste("annot.gene:", annot.gene))
	write.msg(paste("rmasker.filter.margin:", rmasker.filter.margin))
	write.msg(paste("verbose:", verbose))

	rid(sample=sample, dir=dir, chr=chr, no.clipped=no.clipped, ref=ref, rasym=rasym, min.ram=min.ram, 
		oneside.ram=oneside.ram, exo=exo, cbam.chr=cbam.chr, jittering=jittering, 
		merged.rmasker=merged.rmasker, merge.family=merge.family, rmasker.filter.margin=rmasker.filter.margin, 
		annot.oi=annot.oi, annot.gene=annot.gene, verbose=verbose)

} else if (args[1] == "rid.te.primate") {

	chr=NULL; no.clipped=F; rasym=NULL; min.ram=3; oneside.ram=F; exo=F; cbam.chr=T; jittering=2
	merged.rmasker=F; merge.family=F; rmasker.filter.margin=500; annot.oi=T; annot.gene=F; annot.gap=T; verbose=F

	for(i in 2:length(args)) eval(parse(text=args[[i]]))

    write.msg(paste("dir:", dir))
    write.msg(paste("sample:", sample))
    write.msg(paste("chr:", chr))
    write.msg(paste("no.clipped:", no.clipped))
    write.msg(paste("ref:", ref))
    write.msg(paste("rasym:", rasym))
    write.msg(paste("min.ram:", min.ram))
    write.msg(paste("oneside.ram:", oneside.ram))
    write.msg(paste("exo:", exo))
    write.msg(paste("cbam.chr:", F))
    write.msg(paste("jittering:", jittering))
    write.msg(paste("merged.rmasker:", merged.rmasker))
    write.msg(paste("merge.family:", merge.family))
    write.msg(paste("annot.oi:", annot.oi))
    write.msg(paste("annot.gene:", annot.gene))
    write.msg(paste("rmasker.filter.margin:", rmasker.filter.margin))
    write.msg(paste("verbose:", verbose))

    rid.te.primate(sample=sample, dir=dir, chr=chr, no.clipped=no.clipped, ref=ref, rasym=rasym, min.ram=min.ram,
        oneside.ram=oneside.ram, exo=exo, cbam.chr=cbam.chr, jittering=jittering,
        merged.rmasker=merged.rmasker, merge.family=merge.family, rmasker.filter.margin=rmasker.filter.margin,
        annot.oi=annot.oi, annot.gene=annot.gene, verbose=verbose)

} else if (args[1] == "rid.te") { 

	min.ram=3; jittering=2; verbose=F; chr=NULL; no.clipped=F; ref="hg18"; cbam.chr=F

	for(i in 2:length(args)) eval(parse(text=args[[i]]))
	write.msg(paste("sample:", sample))
	write.msg(paste("dir:", dir))
	write.msg(paste("min.ram:", min.ram))
	write.msg(paste("jittering:", jittering))
	write.msg(paste("verbose:", verbose))+   write.msg(paste("chr:", chr))
	write.msg(paste("no.clipped:", no.clipped))
	write.msg(paste("ref:", ref))
	write.msg(paste("cbam.chr:", cbam.chr))

	rid.disc(sample, dir, chr, no.clipped, ref, min.ram, jittering, cbam.chr, verbose)

} else if (args[1] == "somatic") {

	rasym = "um"; matched.control = NULL;
	nonmatched.controls = NULL; ref = "hg18"; oneside.ram = F;
	min.ram = 3; min.acr = 2; 
	min.acrr = 0.5; min.tsd = -15; max.tsd = 30;
	min.score = 0.6; verbose = F; gene.annot = T; ram.cutoff = 6;
	matched.cram = 1; nonmatched.cram = 2; matched.cacr = 1;
	nonmatched.cacr = 2; min.out.conf = 5; no.oi = T; mark.exo = T;
	contig = T; seqid = T
	
	for(i in 2:length(args)) eval(parse(text=args[[i]]))

    write.msg(paste("dir:", dir))
    write.msg(paste("sample:", sample))
	write.msg(paste("ref:", ref))
	write.msg(paste("rasym:", rasym))
	write.msg(paste("min.ram:", min.ram))
	write.msg(paste("min.acr:", min.acr))
	write.msg(paste("oneside.ram:", oneside.ram))
	write.msg(paste("min.acrr:", min.acrr))
	write.msg(paste("min.tsd:", min.tsd))
	write.msg(paste("max.tsd:", max.tsd))
	write.msg(paste("no.oi:", no.oi))
	write.msg(paste("min.out.conf:", min.out.conf))
	write.msg(paste("ram.cutoff:", ram.cutoff))
	write.msg(paste("verbose:", verbose))
	write.msg(paste("contig:", contig))

	call.somatic(dir=dir, sample=sample, ref=ref, rasym=rasym, min.ram=min.ram,
		min.acr=min.acr, oneside.ram=oneside.ram, 
		min.acrr=min.acrr, min.tsd=min.tsd, max.tsd=max.tsd, no.oi=no.oi, min.out.conf=min.out.conf, ram.cutoff=ram.cutoff, 
		verbose=verbose, contig=contig)                     

} else if (args[1] == "germline" || args[1] == "te.primate.germline") {

	ref="hg18"; rasym="ra"; 
	min.ram=3; oneside.ram=F; min.acr=2; min.acrr=0.4; 
	min.tsd=-20; max.tsd=50; no.oi=T; min.out.conf=5; ram.cutoff=6; verbose=F; contig=F

	if (args[1] == "te.primate.germline") {
		no.oi=F; ram.cutoff=3
	}

    for(i in 2:length(args)) eval(parse(text=args[[i]]))

    write.msg(paste("dir:", dir))
    write.msg(paste("sample:", sample))
	write.msg(paste("ref:", ref))
	write.msg(paste("rasym:", rasym))
	write.msg(paste("min.ram:", min.ram))
	write.msg(paste("min.acr:", min.acr))
	write.msg(paste("oneside.ram:", oneside.ram))
	write.msg(paste("min.acrr:", min.acrr))
	write.msg(paste("min.tsd:", min.tsd))
	write.msg(paste("max.tsd:", max.tsd))
	write.msg(paste("no.oi:", no.oi))
	write.msg(paste("min.out.conf:", min.out.conf))
	write.msg(paste("ram.cutoff:", ram.cutoff))
	write.msg(paste("verbose:", verbose))
	write.msg(paste("contig:", contig))

	call.germline(dir=dir, sample=sample, ref=ref, rasym=rasym, min.ram=min.ram,
		min.acr=min.acr, oneside.ram=oneside.ram, 
		min.acrr=min.acrr, min.tsd=min.tsd, max.tsd=max.tsd, no.oi=no.oi, min.out.conf=min.out.conf, ram.cutoff=ram.cutoff, 
		verbose=verbose, contig=contig)
}
