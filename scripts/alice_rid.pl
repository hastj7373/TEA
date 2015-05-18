#!/usr/bin/perl
#
# Alice Eunjung Lee : ejalice.lee@gmail.com
# need to add cbam generation on pq in the pp
# ra.pl cbam /data2/lxyang/cancer_genome/bam_bwa_hg18/ov1411_normal/ov1411_normal.sorted.bam
# final: disc.bam, cl.bam, ram.bam (merged, sorted, index), cbam (index), *.isize, cpos.bz2 

use strict;
use Getopt::Std;
use File::Basename;
use IO::File;

my $usage = qq{
Usage:    rid.pl [options] <sample id>
  
Options: 
         -d STR base_dir (/files/CBMI/parklab/alee/ra)
         -S short reads (no clipped reads and cl.sorted.bam)
         -g STR the whole bam to extract clipped reads 
         -i STR the clipped bam folder
                (e.g. folder/ov0890_cancer.sorted.softclips.consd.sorted.bam, sorted.bam.bai, consd.cpos.bz2)
         -l INT bwa l (40)
         -k INT bwa k (2)
         -n INT bwa n (3)
         -r STR repeat assembly (/files/CBMI/parklab/alee/ra/data/assembly/repeat.combined.div30.isize150.fa)
         -s STR scratch (/scratch/el114)		
         -p STR pq node (ant16.sum)		
         -e execution (NO)		
         -t disable the time stamp (YES)
         -c STR a comma separated chromosome list or species [hs|orangutan] (hs)
         -f force run and rewrite
         -a INT	run the given step only
               	1: rcopying dics.bam and generating indexes and prepare isize and read length (1)
                2: generate ram.bam and ram.bz2 
								3: cbam prep 
                4: run rid per chr and merged them
         -b STR rcopying *disc.bam, *cl.disc.bam, *isinfo from the given dir
                (e.g.) -b alee\@ant13.sum:~.... then rcopy files from ant13.sum to the base_dir
         -m LSF_based parallel runs for each chr (NO)

};

my %opts = ();
getopts("hSc:s:q:d:l:k:n:r:s:p:a:etb:fg:i:m", \%opts);

if ( @ARGV < 1 || defined($opts{h}) ) { die $usage }
my $sample = shift;

# set parameter default values
my $q = "all_unlimited";
my $bdir = "/files/CBMI/parklab/alee/ra2";
my $max_mem = "2000000000";
my $l = 40;
my $k = 2;
my $n = 3;
my $ra = "/files/CBMI/parklab/alee/ra/data/assembly/repeat.combined.div30.isize150.fa";
my $sourcedir = "";
#my $sourcedir = "alee\@ant16.sum:~lxyang/cancer_genome/sv_meerkat/run";
my $scratch = "/scratch/el114";
my $pq = "ant16.sum";
my $exec = 0;
my $time = 1;
my $chr = "hs";
my $step=  -1;
my $force =  0;
my $bamf = "";
my $cbamdir = "";
my $short = 0;
my $parallel = 0;

$sourcedir = $opts{b} if (defined($opts{b}));
$bdir = $opts{d} if (defined($opts{d}));
$l = $opts{l} if (defined($opts{l}));
$k = $opts{k} if (defined($opts{k}));
$n = $opts{n} if (defined($opts{n}));
$ra = $opts{r} if (defined($opts{r}));
$scratch = $opts{s} if (defined($opts{s}));
$pq = $opts{p} if (defined($opts{p}));
$chr = $opts{c} if (defined($opts{c}));
$step = $opts{a} if (defined($opts{a}));
$bamf = $opts{g} if (defined($opts{g}));
$cbamdir = $opts{i} if (defined($opts{i}));
$exec = 1 if (defined($opts{e}));
$time = 0 if (defined($opts{t}));
$force = 1 if (defined($opts{f}));
$short = 1 if (defined($opts{S}));
$parallel = 1 if (defined($opts{m}));

my ($dir, $ts, $time, $cmd) = ("$bdir/$sample", "", "", "");
if ($time) { $ts = "\\time -f \"%E elapsed\""; } 

my $time_start = time;
###############################################
# step1: copying discordant bams and isize files from sourcedir to dir 
if (! -d "$dir/bam") { 		
	$cmd = "mkdir -p $dir/bam";
	print $cmd."\n";
	if ($exec) { system("$cmd") }
	print "echo \"$dir/bam created\"\n";
}

if ($sourcedir ne "" & ($step == -1 || $step == 1)) {
	if (! -e "$dir/bam/$sample.disc.bam" && !-e "$dir/bam/$sample.disc.sorted.bam" || $force) {
		print "echo \"(r)copying disc.sorted.bam...\"\n";
		if ($sourcedir =~ m/\@/) {
			$cmd = "$ts rcp $sourcedir/$sample.disc.sorted.bam $dir/bam";
			print $cmd."\n";
			if ($exec) { system("$cmd") }
		} else {
			$cmd = "$ts cp $sourcedir/$sample.disc.sorted.bam $dir/bam";
			print $cmd."\n";
			if ($exec) { system("$cmd") }
		}
		print "echo \"done (r)copying disc.sorted.bam.\"\n";
	}

	if (!$short && ! -e "$dir/bam/$sample.cl.disc.bam" && !-e "$dir/bam/$sample.cl.disc.sorted.bam" || $force) {
		print "echo \"(r)copying cl.disc.sorted.bam...\"\n";
		if ($sourcedir =~ m/\@/) {
			$cmd = "$ts rcp $sourcedir/$sample.cl.disc.sorted.bam $dir/bam";
			print $cmd."\n";
			if ($exec) { system("$cmd") }
		} else {
			$cmd = "$ts cp $sourcedir/$sample.cl.disc.sorted.bam $dir/bam";
			print $cmd."\n";
			if ($exec) { system("$cmd") }
		}
		print "echo \"done (r)copying cl.disc.sorted.bam.\"\n";
	}

	if (! -e "$dir/bam/$sample.isinfo" || $force) {
		print "echo \"(r)copying isinfo...\"\n";
		if ($sourcedir =~ m/\@/) {
			$cmd = "$ts rcp $sourcedir/$sample.isinfo $dir/bam";
			print $cmd."\n";
			if ($exec) { system("$cmd") }
		} else {
			$cmd = "$ts cp $sourcedir/$sample.isinfo $dir/bam";
			print $cmd."\n";
			if ($exec) { system("$cmd") }
		}
		print "echo \"done (r)copying isinfo.\"\n";
	}
} elsif ($sourcedir eq "") {
	print "echo \"skipping copying disc.bams: sourcedir was not specified.\"\n";
}

if (! -e "$dir/bam/$sample.disc.bam.bai") {
	print "echo \"generating the indexe for disc.bam...\"\n";
	$cmd = "samtools sort $dir/bam/$sample.disc.sorted.bam $dir/bam/$sample.disc; samtools index $dir/bam/$sample.disc.bam";
	print $cmd."\n";
	if ($exec) { system("$cmd") }
	print "echo \"done generating the index for $dir/bam/$sample.disc.bam.\"\n";
}
if (!$short && ! -e "$dir/bam/$sample.cl.disc.bam.bai") {
	print "echo \"generating the indexe for cl.disc.bam...\"\n";
	$cmd = "samtools sort $dir/bam/$sample.cl.disc.sorted.bam $dir/bam/$sample.cl.disc; samtools index $dir/bam/$sample.cl.disc.bam";
	print $cmd."\n";
	if ($exec) { system("$cmd") }
	print "echo \"done generating index for $dir/bam/$sample.cl.disc.bam.\"\n";
}

if (! -e "$dir/bam/$sample.isize") {
 	print "echo \"formatting isize and rl information\"\n";
  $cmd = "ra.pl cisize $dir/bam/$sample.isinfo";
  print $cmd."\n";
  if ($exec) { system($cmd) }
 	print "echo \"done formatting isize and rl.\"\n";
}

if (-e "$dir/bam/$sample.disc.sorted.bam" || -e "$dir/bam/$sample.cl.disc.sorted.bam" || -e "$dir/bam/$sample.isinfo") {
	$cmd = "rm $dir/bam/$sample.disc.sorted.bam $dir/bam/$sample.cl.disc.sorted.bam $dir/bam/$sample.isinfo";
	print $cmd."\n";
	if ($exec) { system($cmd) }
	print "echo \"done rm disc.sorted.bam, cl.disc.sorted.bam, isinfo\"\n";
}

if ($step == 1) { exit }

###############################################
# step2: generate the merged ram.bam and ram.bz2 
if ($step == -1 || $step == 2) {
	if (! -e "$dir/bam/$sample.ram.bam" || $force) {
		if (! -e "$dir/bam/$sample.disc_1.fastq.gz" || !-e "$dir/bam/$sample.disc_2.fastq.gz") {
			print "echo \"generating fastqs from disc.bam...\"\n";
			$cmd = "$ts ra.pl bam2fq -p $dir/bam/$sample.disc.bam $dir/bam";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating fastqs from disc.bam.\"\n";
		}
		if (!$short && (!-e "$dir/bam/$sample.cl.disc_1.fastq.gz" || !-e "$dir/bam/$sample.cl.disc_2.fastq.gz")) {
			print "echo \"generating fastqs from cl.disc.bam...\"\n";
			$cmd = "$ts ra.pl bam2fq -p $dir/bam/$sample.cl.disc.bam $dir/bam";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating fastqs from cl.disc.bam.\"\n";
		}

		my @fql = ();
		push(@fql, "$sample.disc_1.fastq.gz");
		push(@fql, "$sample.disc_2.fastq.gz");
		if (!$short) {
			push(@fql, "$sample.cl.disc_1.fastq.gz");
			push(@fql, "$sample.cl.disc_2.fastq.gz");
		}

		my $sra = "/scratch/el114/assembly/".basename($ra);
		if (!-e "/scratch/el114/assembly") {
			system("mkdir -p /scratch/el114/assembly")
		}
		if (!-e $sra) {
			$cmd = "cp /files/CBMI/parklab/alee/ra/data/assembly/".basename($ra)."* /scratch/el114/assembly/";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"copying assembly files to scratch\"\n";
		}
		for my $fq (@fql) {
			my $prefix = $fq;
			$prefix =~ s/.fastq.gz//;

			print "echo \"aligning $fq to $sra...\"\n";
			$cmd = "$ts bwa aln -l $l -k $k -n $n $sra $dir/bam/$fq > $scratch/$prefix.sai; bwa samse -n 100 $sra $scratch/$prefix.sai $dir/bam/$fq | samtools view -bt $sra.fai - > $dir/bam/$prefix.ra.bam";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done alignment for $fq.\"\n";
		}

		if (!-e "$dir/bam/$sample.disc.ram.bam") {
			print "echo \"generating ram.bam for disc.bam\"\n";
			$cmd = "ra.pl bamram -o $dir/bam $dir/bam/$sample.disc.bam $dir/bam/$sample.disc_1.ra.bam $dir/bam/$sample.disc_2.ra.bam";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating ram.bam for disc.bam\"\n";
		} else {
			print "echo \"$dir/bam/$sample.disc.ram.bam already exists\"\n";
		}

		if (!$short && !-e "$dir/bam/$sample.cl.disc.ram.bam") {
			print "echo \"generating ram.bam for cl.disc.bam\"\n";
			$cmd = "ra.pl bamram -o $dir/bam $dir/bam/$sample.cl.disc.bam $dir/bam/$sample.cl.disc_1.ra.bam $dir/bam/$sample.cl.disc_2.ra.bam"; 
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating ram.bam for cl.disc.bam\"\n";
		} else {
			print "echo \"$dir/bam/$sample.cl.disc.ram.bam already exists\"\n";
		}
			
		if ($short) {
			$cmd = "mv $dir/bam/$sample.disc.ram.bam  $dir/bam/$sample.ram.raw.bam; samtools sort $dir/bam/$sample.ram.raw.bam $dir/bam/$sample.ram; samtools index $dir/bam/$sample.ram.bam; rm $dir/bam/$sample.ram.raw.bam";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating sorted ram.bam.\"\n";

			$cmd = "mv $dir/bam/$sample.disc.ram.bz2 $dir/bam/$sample.ram.bz2";
			print $cmd."\n";
			if ($exec) { system($cmd) }
		} else {
			print "echo \"merging disc.ram.bam cl.disc.ram.bam...\"\n";
			$cmd = "samtools merge $dir/bam/$sample.ram.raw.bam $dir/bam/$sample.disc.ram.bam $dir/bam/$sample.cl.disc.ram.bam; samtools sort $dir/bam/$sample.ram.raw.bam $dir/bam/$sample.ram; samtools index $dir/bam/$sample.ram.bam; rm $dir/bam/$sample.ram.raw.bam";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating the merged ram.bam.\"\n";

			print "echo \"merging disc.ram.bz2 and cl.disc.ram.bz2...\"\n";
			$cmd = "bzcat $dir/bam/$sample.disc.ram.bz2 > $dir/bam/$sample.ram; bzcat $dir/bam/$sample.cl.disc.ram.bz2 >> $dir/bam/$sample.ram; bzip2 $dir/bam/$sample.ram";
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done generating the merged.ram.bz2.\"\n";
		}
	} else {
		print "echo \"$dir/bam/$sample.ram.bam already exists.\"\n";
	}
}

# clean-up
#$cmd = " rm $dir/bam/*.fastq.gz $dir/bam/*.ra.bam";
#print $cmd."\n";
#if ($exec) { system($cmd) }

if ($step == 2) { exit }
 
###############################################
# step3: cbam prep
if (!$short && ($step == -1 || $step == 3)) {
	if ($cbamdir ne "") {
		$cmd = "cp $cbamdir/$sample.sorted.softclips.consd.* $dir/bam";
    print $cmd."\n";
    if ($exec) { system($cmd) }
		print "echo \"done copying the clipped bam from $cbamdir\"\n";
	} elsif ($bamf ne "") {
		$cmd = "ra.pl cbam -o $dir/bam $bamf";
    print $cmd."\n";
    if ($exec) { system($cmd) }
		print "echo \"done generating the clipped bam from $bamf\"\n";
	} elsif (! -e "$dir/bam/$sample.sorted.softclips.consd.sorted.bam" &&  !-e "$dir/bam/$sample.sorted.softclips.consd.bam") {
			die "cbam files in $dir/bam, [-i cbamdir], or [-g bamf] are required";
	} elsif (-e "$dir/bam/$sample.sorted.softclips.consd.sorted.bam" ||  -e "$dir/bam/$sample.sorted.softclips.consd.bam") {
		print "echo \"the clipped file already exists in $dir/bam\"\n";
	}
}
if ($step == 3) { exit }

###############################################
# step4: run rid per chr
if ( $step == -1 || $step == 4 ) {
	#if ($chr eq "all") {
	#	print "#running rid for all chromosomes...\n";
	#	$cmd = "$ts R --vanilla --args rid.disc sample=\\\"$sample\\\" < ~/repeat_analysis/code/run.rid.r";
	#	print $cmd."\n";
	#	if ($exec) { system($cmd) }
	#	print "#done rid.\n\n";
	#} else {
	my @chrl = ();
	if ($chr eq "hs") { 
		for (my $i=1; $i<=22; $i++) { push(@chrl, "chr$i"); } push(@chrl, "chrX"); push(@chrl, "chrY"); 
	} elsif ($chr eq "orangutan") {
		for (my $i=1; $i<=22; $i++) { 
			if ($i==2) { push(@chrl, "chr2a"); push(@chrl, "chr2b") }	
			else { push(@chrl, "chr$i"); } 
		}
		push(@chrl, "chrX");
	} elsif ($chr =~ m/,/) {	
		@chrl = split(/,/, $chr);
	}
		for my $chr (@chrl) {
			if (! -e "$dir/cluster/$sample.$chr.cluster") {
			print "echo \"running rid for $chr ...\"\n";
			if ($short) {
				if ($parallel) {
					$cmd = "bsub -q $q -J $sample.$chr -o log/$sample.$chr.log \"R --vanilla --args rid.disc dir=\\\\\\\"$bdir\\\\\\\" sample=\\\\\\\"$sample\\\\\\\" chr=\\\\\\\"$chr\\\\\\\" short=T < ~/repeat_analysis/code/run.rid.r\"";
				} else {
					$cmd = "$ts R --vanilla --args rid.disc dir=\\\"$bdir\\\" sample=\\\"$sample\\\" chr=\\\"$chr\\\" short=T < ~/repeat_analysis/code/run.rid.r";
				}
			} else {
				if ($parallel) {
					$cmd = "bsub -q $q -J $sample.$chr -o log/$sample.$chr.log \"R --vanilla --args rid.disc dir=\\\\\\\"$bdir\\\\\\\" sample=\\\\\\\"$sample\\\\\\\" chr=\\\\\\\"$chr\\\\\\\" < ~/repeat_analysis/code/run.rid.r\"";
				} else {
					$cmd = "$ts R --vanilla --args rid.disc dir=\\\"$bdir\\\" sample=\\\"$sample\\\" chr=\\\"$chr\\\" < ~/repeat_analysis/code/run.rid.r";
				}
			}
			print $cmd."\n";
			if ($exec) { system($cmd) }
			print "echo \"done rid for $sample, $chr.\"\n";
			}
		}
		my $done = 1;
		for my $chr (@chrl) {
			if (! -e "$dir/cluster/$sample.$chr.cluster") { $done = 0; last }
		}
		if ($done) { 
			$cmd = "echo `date` > $dir/cluster/rid.ok"; 
			print $cmd . "\n";	
			if ($exec) { system($cmd) }
		}
}

# running time stat
if ($exec) {
	my $time_end = time;
	my $diff = $time_end - $time_start;
	my $sec  =  $diff % 60;
	$diff = ($diff - $sec) / 60;
	my $min =  $diff % 60;
	$diff = ($diff - $min) / 60;
	print "Total Elapsed Time: $diff:$min:$sec\n\n"; 
}
