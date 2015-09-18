#!/usr/bin/perl

use strict;
no strict "refs";
use Getopt::Std;
use File::Basename;
use IO::File;
use IO::Compress::Gzip qw(gzip $GzipError);
use List::Util qw(min max sum);
#use warnings;
my $cmd = basename($0);

# user-defined variables
#my $blast_path = "/home/el114/download/ncbi-blast-2.2.26+/bin";
#my $blast_path = "/bitools/BLAST/ncbi-blast-2.2.27+/bin";
#my $bwa = "/bitools/BWA/bwa-0.7.5a/bwa";
#my $samtools = "/bitools/SAMTools/samtools-0.1.19/samtools";
#my $assembler = "/bitools/CAP3/CAP3_07-10-15/cap3";

my $line;
my %hash = ();

open my $f, "<", "conf_file";

while($line=<$f>){
     my @val = split('=', $line);

     $val[0] =~ s/^\s+//;
     $val[0] =~ s/\s+$//;
     $val[1] =~ s/^\s+//;
     $val[1] =~ s/\s+$//;
     $hash{$val[0]} = $val[1];
}
close $f;

my $blast_path = "$hash{'blast_home'}";
my $bwa = "$hash{'bwa_home'}"; $bwa = "$bwa/bwa";
my $samtools = "$hash{'samtools_home'}"; $samtools = "$samtools/samtools";
my $assembler = "$hash{'cap3_home'}"; $assembler = "$assembler/cap3";
my $bzip2 = "bzip2";
my $bedtools = "$hash{'bedtools_home'}"; $bedtools = "$bedtools/bedtools";
my $reference = "$hash{'fasta_location'}";
my $retro_ref = "$hash{'retro'}";
my $temp_loc = "$hash{'tmp_dir'}";


my $usage = qq{

  Usage:   $cmd <command> [options]

  Command: tea    generate a shell script to run tea
  
           bam2fq  generate fastq files from a bam file
  
           cbam	   extract the soft-clipping positions and make a bam only with the clipped reads

           bamram  extract the coordinates of rams  and generate a bam only with rams

           isize   extract the distribution of fragment sizes and write it to the standard output

           cisize  format the isize and read length information from meerkat output
            
           contig  generate contigs of ram mates and clipped sequences

           map_consensus  <infile> <family:L1|Alu|PABL_A-INT> <map_subfmaily: 1|0> 

           rmasker_coordinates <input repeat masker file: *.gz> <1st record line#> <verbose: 1|0>
           rmasker_coordinates_strand <input repeat masker file: *.gz> <1st record line#> <verbose: 1|0>

           ram     extract the coordinates of rams per te type

           sram    separte rams per each repeat type

           mapcnt  count the mapped reads for each sequence
 
           rmdup    remove PCR duplicates in a bam file with unmapped reads

};

my $log_dir = ".";

if (@ARGV < 1) { die $usage }

$cmd = shift;
my @functions = qw(tea bam2fq cbam ram bamram isize cisize sram contig map_consensus make_subfamilyseq 
                   rmasker_coordinates mapcnt rmdup rmasker_coordinates_strand);

if ($cmd eq "map_consensus") {
	my ($infile, $family, $map_subfamily) = @ARGV;
	&map_consensus($infile, $family, $map_subfamily, 0);
	exit(1);
} 

if ($cmd eq "make_subfamilyseq") {
	my ($family) = @ARGV;
	&make_subfamilyseq($family);
	exit(1);
}

if (grep {$_ eq $cmd} @functions) { &$cmd(); }
else { die $usage; }

sub tea 
{
  my $usage = qq{
  Usage:    tea.pl [options] <sample id>
            generate a shell script to perform the following steps 
            1. procss input files--sorting, indexing, an formatting
               disc.sorted.bam 
               cl.disc.sorted.bam (if no.clipped=F)
               cl.sorted.bam (if exogenous=T)
               isinfo (isize and read length) 
             2. generate fastq files
             3. mapping
             4. generate ram.bam and ram.bz2
             5. cbam prep
             6. run (rid | rid.te | rid.te.primate) 
             7. run (germline | te.primate.germline)
	  
  Options: -d STR  base dir (.)
           -p INT  number of threads for bwa alignment          
           -P STR  tea base dir to add to \$PATH
           -F      skip the step four, cbam generation (no)
           -g STR  a bam file to extract clipped reads
           -l INT  bwa aln seed lenegh (40)
           -k INT  bwa aln maximum differences in the seed (2)
           -n INT  bwa aln max #diff (int) or missing prob under 0.02 err rate (floot) (3)
           -r STR  sequence assembly (.fasta)
           -f      family-wise merge of rams (no)
           -R STR  assembly symbol (ra and va for endogenoeus and exogenous)
           -s STR  scratch or temp space to save the intermediate mapping files (/tmp)
           -c STR  a comma separated chromosome list or ref
                   [hg18|hg19|panTro3|ponAbe2|rheMac2] (hg18)
           -a INT  start from the given step
           -x      exogenous (no)
           -S      no clipped reads and cl.sorted.bam
           -C      cleaning up the intermediate files such as fastq, ra.bam, etc. (no)
           -M INT  minimum ram
           -o      also reporting clusters with only one side ram cluster (no)
           -j INT  jittering or allowed gap in extracting clipped positions (2)
           -e STR  rid execution function (rid | rid.te | rid.te.primate) (rid)
           -t      contig generation (no) 
		   -T	transduction (no)
		   -O	orphan (no)
           ### LSF parameters 
           -m      LSF job for each chr (no)
           -L STR  LSF log dir (./log)
           -q STR  LSF queue (long)
           -b STR  LSF estimated run time (72:00)

	};

	my %opts = ();
	getopts("hd:g:l:k:n:r:s:c:Smq:GP:a:p:xR:j:oOM:fFie:tTb:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage };

	my $sample = shift(@ARGV);
	
	# set defaults
	my $bdir = ".";
	my $bamf = "";
	my $logdir = "./log";
	my $l = 40;
	my $k = 2;
	my $n = 3;
	my $ref = "hg18";
	my $no_clipped = 0;
	my $parallel = 0;
	my $tea = ""; 
	my $scratch = "/tmp";
	my $q = "long";
	my $sortmem = "2000000000";
	my $clean = 0;
	my $step = 1;
	my $nproc = 1;
	my $exogenous = 0;
	my $minram = 3;
	my $oneside_ram = 0;
	my $jittering = 2;
	my $rasym = "ra";
	my $family = 0;
	my $skip_cbam = 0;
	my $exec = "rid";
	my $contig = 1;
	my $time = "72:00";
	my $orphan = 0;
	my $transduction = 0;
	my $merges = 0;


	$bdir = $opts{d} if (defined($opts{d}));
	my $dir = "$bdir/$sample";
	$bamf = $opts{g} if (defined($opts{g}));
	$l = $opts{l} if (defined($opts{l}));
	$k = $opts{k} if (defined($opts{k}));
	$n = $opts{n} if (defined($opts{n}));
	$ref = $opts{c} if (defined($opts{c}));
	$exec = $opts{e} if (defined($opts{e}));
	$no_clipped = 1 if (defined($opts{S}));
	$parallel = 1 if (defined($opts{m}));
	$clean = 1 if (defined($opts{C}));
	$q = $opts{q} if (defined($opts{q}));
	$logdir = $opts{L} if (defined($opts{L}));
	$step = $opts{a} if (defined($opts{a}));
	$nproc = $opts{p} if (defined($opts{p}));
	$exogenous = 1 if (defined($opts{x}));
	$minram = $opts{M} if (defined  $opts{M});
	$oneside_ram = 1 if (defined $opts{o});
	$family = 1 if (defined $opts{f});
	$jittering = $opts{j} if (defined $opts{j});
	$skip_cbam = 1 if (defined $opts{F});
	$contig = 0 if (defined $opts{t});
	$time = $opts{b} if (defined $opts{b});

	$orphan = 1 if (defined $opts{O});
	$transduction = 1 if (defined $opts{T});
	##HC
	$merges = 1 if (defined $opts{G});

	my $ra = "$tea/lib/assembly/repeat.combined.div30.isize150.fa";	
  
	$ra = $opts{r} if (defined($opts{r}));
	$rasym = $opts{R} if (defined($opts{R}));
	$scratch = $opts{s} if (defined($opts{s}));
	$tea = $opts{P} if (defined($opts{P}));
  #my $ra = "$tea/lib/assebmly/repeat.combined.div30.isize150.fa";
  my $ra_um = "$tea/lib/assembly/repeat.combined.div30.isize150.fa";
	my $bash_binary = `which bash`;
	print "#!$bash_binary\n";
	print "export tea_base=$tea\n";
	
    #Only harvard	
 	print "export PATH=$tea/R:$tea/scripts:\$PATH\n"; 
	print "module load stats/R/2.15.3\n";
  print "export R_LIBS=/opt/R-2.15.3/lib64/R/library/:/home/el114/hyunchul/R/library\n";
	print "export PERL5LIB=/home/ly55/perl/share/perl/5.10.0:/home/ly55/perl/lib/perl5/:/home/ly55/perl/lib/perl/5.10.0:/home/el114/hyunchul/lib/perl5/:/home/el114/hyunchul/lib/share/perl5/\n";




	my $ascratch = "$scratch/assembly";
	my $sscratch = "$scratch/$sample";
	
	################################################################################
	#  Note that samtools and tea.pl must be in $PATH for the code below to work.  
	#  Tea requires the below input files in $dir/bam.  
	#    sample.disc.sorted.bam (if exogenous=F)
	#    sample.cl.disc.sorted.bam (if exogenous=F & no_clipped=F)
	#    sample.cl.sorted.bam (if exoengeous=T & no_clipped=F)
	#    sample.mapped_um.sorted.bam (if exoengeous=T & no_clipped=T)
	#    sample.isinfo

	my $bamdir = "$dir/bam";
	my $cldir = "$dir/cluster\_${rasym}m";
	my $dbamf = "$bamdir/$sample.disc.sorted.bam";
	my $dbamf2 = "$bamdir/$sample.disc.bam";
	my $clbamf = "$bamdir/$sample.cl.sorted.bam"; # a bam with unmapped and soft-cilpped reads 
	my $umbamf = "$bamdir/$sample.mapped_um.sorted.bam"; # a bam with unmapped reads
	my $umbamf1 = "$bamdir/$sample.um.raw.bam";
	my $umbamf2 = "$bamdir/$sample.um.bam";
	my $umbz2 = "$bamdir/$sample.umm.bz2";
	my $cldbamf = "$bamdir/$sample.cl.disc.sorted.bam"; 
	my $cldbamf2 = "$bamdir/$sample.cl.disc.bam"; 
	my $isf = "$bamdir/$sample.isinfo";
	my $isizef = "$bamdir/$sample.isize";
	my $rlf = "$bamdir/$sample.rl";
	my $cbamf = "$bamdir/$sample.softclips.consd.bam";
	my $cbamf2 = "$bamdir/$sample.sorted.softclips.consd.bam";
	my $fprefix = "";

	if ($step <= 1) {
		if ($exogenous && $no_clipped) {
			if (!-e $umbamf) { 
				print STDERR "no bam with unmapped reads: $umbamf";
				die "no bam with unmapped reads: $umbamf";
			}
			else { 
				print "ln -s $umbamf $umbamf2; ln -s $umbamf.bai $umbamf2.bai\n"; 
			}
		} elsif ($exogenous && !$no_clipped) {
			if (!-e $clbamf) { 
				print STDERR "no bam with unmapped and soft clipped reads: $clbamf"; 
				die "no bam with unmapped and soft clipped reads: $clbamf"; 
			} #elsif (-e $umbamf2) {
				#print STDERR "$umbamf2 already exists";
				#die "$umbamf2 already exists";} 
				else {
				# better make a function for the task below: generating um.bam and umm.bz2
				
				if ($transduction)  {
					print qq{
						python $tea/scripts/python/transduction_module.py $bdir $sample $tea $bedtools $samtools
					}
				}
				if ($orphan )  {
					$contig=1;
					print qq{
						python $tea/scripts/python/orphan_module.py $bdir $sample $tea $bedtools $samtools
					}
				}
				if ($transduction == 0 && $orphan == 0){
				print qq{
				echo "generating a bam with unmapped reads ..."
				$samtools view -hX $clbamf | perl -ne '\@a=split(/\\t/); if (m/^@/ || \$a[1] =~ /[uU]/) { print } ' | $samtools view -bS - > $umbamf1 
				tea.pl rmdup $umbamf1 $umbamf2
				$samtools index $umbamf2
				echo "generating unmapped read positions umm.bz2 ..."
				#$samtools view -X $umbamf2 | perl -e 'my %h = (); while (<>) { \@a=split(/\\t/); if (\$a[1] =~ m/U/ && !(\$a[1] =~ m/u/) && \$a[4]>0 && (\$_ =~ /XT:A:U/ || !(\$_ =~ /XT:A/))) { \$r = \$a[0]; \$r =~ s/sc\$//; \$r =~ s/mu[12]\$//;  if (!exists(\$h{\$r})) { \$a[2] =~ s/chr//; \$h{\$r} = 1; if (\$a[1] =~ /r/) { print "\$a[0]\\t\$a[2]\\t-\$a[3]\\tx\\n"} else { print "\$a[0]\\t\$a[2]\\t\$a[3]\\tx\\n"} } }}' | $bzip2 - > $umbz2 
				#$samtools view -X um2_sorted.bam | perl -e 'my %h = (); while (<>) { \@a=split(/\\t/); if (\$a[1] =~ m/U/ && !(\$a[1] =~ m/u/) && \$a[4]>0 && (\$_ =~ /XT:A:U/ || !(\$_ =~ /XT:A/))) { \$r = \$a[0]; \$r =~ s/sc\$//; \$r =~ s/mu[12]\$//;  if (!exists(\$h{\$r})) { \$a[2] =~ s/chr//; \$h{\$r} = 1; if (\$a[1] =~ /r/) { print "\$a[0]\\t\$a[2]\\t-\$a[3]\\tx\\n"} else { print "\$a[0]\\t\$a[2]\\t\$a[3]\\tx\\n"} } }}' | $bzip2 - > $umbz2 
				}
				}
	
				}
		} else {
			if (!-e $dbamf) { die "no $dbamf in $bamdir" }
			print "echo \"sorting and generating the index for $dbamf...\"\n";
			$fprefix = $dbamf2;
			$fprefix =~ s/\.bam$//;
			print "$samtools sort -m $sortmem $dbamf $fprefix\n";
      print "$samtools index $fprefix.bam\n";
	
			if (!$no_clipped) {
				if (!-e $cldbamf) { die "no $cldbamf in $bamdir" }
			  print "echo \"sorting and generating the index for $cldbamf...\"\n";
				$fprefix = $cldbamf2;
				$fprefix =~ s/\.bam$//;
			  print "$samtools sort -m $sortmem $cldbamf $fprefix\n";
        print "$samtools index $fprefix.bam\n";
			}
		}
	
		if (!-e $isizef || !-e $rlf) {	
			if (!-e $isf) { die "$isf is required" } 
	  	print "echo \"formatting isize and rl information...\"\n";
			print "tea.pl cisize $isf\n";
		}

		# clean-up
		if ($clean) {
 	 		print "echo \"removing the input files...\"\n";
			print "rm -f $dbamf $clbamf $cldbamf $isf\n";
		}
	}
	
	#  These files should now exist:
	#    sample.disc.bam        (if !exogenous)
	#    sample.disc.bam.bai    (if !exogenous)
	#    sample.cl.disc.bam     (if !exogenous && !no_clipped)
	#    sample.cl.disc.bam.bai (if !exogenous=F && !no_clipped) 
	#    sample.cl.bam          (if exogenous)
	#    sample.cl.bam.bai      (if exogenous)
	#    sample.isize
	#    sample.rl
	################################################################################
	
	if ($step <=2 && $rasym ne "um") {
  	if ($exogenous) {
			print "echo \"generating fastqs from $umbamf2...\"\n";
			print "tea.pl bam2fq -p $umbamf2 $bamdir\n";
		} else {
			print "echo \"generating fastqs from $dbamf2...\"\n";
			print "tea.pl bam2fq -p $dbamf2 $bamdir\n";
			if (!$no_clipped) {
		 		print "echo \"generating fastqs from $cldbamf2...\"\n";
				print "tea.pl bam2fq -p $cldbamf2 $bamdir\n";
			}
		}
	}

	my @fql = ();
	if ($exogenous) {
	  push(@fql, "$sample.um_1.fastq.gz");
	  push(@fql, "$sample.um_2.fastq.gz");
	} else {
	  push(@fql, "$sample.disc_1.fastq.gz");
	  push(@fql, "$sample.disc_2.fastq.gz");
	  if (!$no_clipped) {
	    push(@fql, "$sample.cl.disc_1.fastq.gz");
	    push(@fql, "$sample.cl.disc_2.fastq.gz");
	  }
	}

	if ($step <=3 && $rasym ne "um") {
		if (! -e $ascratch) { print "mkdir -p $ascratch\n"; } # assembly scratch
		if (! -e $sscratch) { print "mkdir -p $sscratch\n"; } # sample scratch

		if (!-e $ascratch."/".basename($ra).".bwt") {
			print "echo \"copying bwa indexes to $ascratch\"\n";
		 	print "cp $ra* $ascratch\n";
		}
		$ra = $ascratch."/".basename($ra);
	
		for my $fq (@fql) {
			my $prefix = $fq;
		  	$prefix =~ s/.fastq.gz//;
			$prefix .= ".$rasym";

			my $ffq = "$bamdir/$fq"; 
			my $sai = "$sscratch/$prefix.sai";
			my $fai = "$ra.fai";
			my $rbam = "$bamdir/$prefix.bam";
		
		 	print qq{
			echo "aligning $fq to $ra..."
			$bwa aln -t $nproc -l $l -k $k -n $n $ra $ffq > $sai
			$bwa samse -n 100 $ra $sai $ffq | $samtools view -bt $fai - > $rbam

			};
		}
	}

	if ($step <=4 && $rasym ne "um") {
	  # generate bams with rams (*.ram.bam) and files with ram positions (*.ram.bz2)\"\n";
		my (@rbaml, @rambaml, @rambz2l, $rbam_prefix, $refbam, $rbam1, $rbam2, $rambam, $rambz2);

		if ($exogenous) { $rbam_prefix ="$sample.um" }
		else { $rbam_prefix ="$sample.disc" }

		$refbam = "$rbam_prefix.bam";
		$rbam1 = "$rbam_prefix\_1.$rasym.bam";
		$rbam2 = "$rbam_prefix\_2.$rasym.bam";
		$rambam = "$rbam_prefix.${rasym}m.bam"; # ram or vam
		$rambz2 = "$rbam_prefix.${rasym}m.bz2";
		print "tea.pl bamram -d $bamdir $refbam $rambam $rambz2 $rbam1 $rbam2\n";
		push(@rbaml, "$bamdir/$rbam1"); # for clean-up
		push(@rbaml, "$bamdir/$rbam2");
		push(@rambaml, "$bamdir/$rambam");
		push(@rambz2l, "$bamdir/$rambz2");

		if (!$exogenous && !$no_clipped) {
			$rbam_prefix ="$sample.cl.disc";
			$refbam = "$rbam_prefix.bam";
			$rbam1 = "$rbam_prefix\_1.$rasym.bam";
			$rbam2 = "$rbam_prefix\_2.$rasym.bam";
			$rambam = "$rbam_prefix.${rasym}m.bam";
			$rambz2 = "$rbam_prefix.${rasym}m.bz2";
			print "tea.pl bamram -d $bamdir $refbam $rambam $rambz2 $rbam1 $rbam2\n";
			push(@rbaml, "$bamdir/$rbam1"); # for clean-up
			push(@rbaml, "$bamdir/$rbam2");
			push(@rambaml, "$bamdir/$rambam");
			push(@rambz2l, "$bamdir/$rambz2");
		}
	
		# merge ram.bam files and ram.bz2 files	
		my $temp = "$bamdir/$sample.${rasym}m.raw.bam";
		my $m_rambam_prefix = "$bamdir/$sample.${rasym}m";
		my $m_rambam = "$m_rambam_prefix.bam";

		if (@rambaml >1) {
			print "\t\t$samtools merge $temp @rambaml\n";
		} else {
			print "\t\tmv $rambaml[0] $temp\n";
		}
		print qq{
			$samtools sort -m $sortmem $temp $m_rambam_prefix
			$samtools index $m_rambam
			rm $temp
		};

		if (@rambz2l >1) {
			print qq{
				rm -f $m_rambam_prefix.bz2 
		 		bzcat @rambz2l >> $m_rambam_prefix
				bzip2 -f $m_rambam_prefix
			};
		} else {
			print "\t\tmv $rambz2l[0] $m_rambam_prefix.bz2\n";
		}

		# clean-up
		if ($clean) { print "rm -f @fql @rbaml @rambaml @rambz2l\n" }
	}
	
	# cbam prep
	if ($step <=5 && !$skip_cbam ) {
 		if ($bamf ne "") {
			print "tea.pl cbam -o $bamdir $bamf\n";
 	 	} elsif (! -e $cbamf) {
		  die "no $cbamf and no input bam specified to make a clipped bam [-g]";
		}
	} elsif ($skip_cbam) {
		print "echo \"skipping cbam generation\"\n";
	}

	# run rid
	if ($step <=6) {

		my @chrl = ();
		my %chrl = ();
		&get_chrl(\@chrl, \%chrl, $ref);

		my $R = "R --no-save --no-restore --no-environ --args $exec";

		my $booleanS = "";
		# test whether cbam has a 'chr' prefix; for compatibility, check a cbam with 'sorted' infix 
		my $cbam_chr = 0;
		if (-e $cbamf) {  $cbam_chr = &bam_chr($cbamf) }
		elsif (-e $cbamf2) { $cbam_chr = &bam_chr($cbamf2) }
		if ($cbam_chr) { 
			$booleanS = "cbam.chr=T" 
		} else {
			$booleanS = "cbam.chr=F" 
		}

		if ($no_clipped) { $booleanS .= " no_clipped=T" } else { $booleanS .= " no_clipped=F" }
		if ($oneside_ram) { $booleanS .= " oneside.ram=T" } else { $booleanS .= " oneside.ram=F" }
		if ($exogenous) { $booleanS .= " exo=T" } else {  $booleanS .= " exo=F" }
		if ($family) { $booleanS .= " merge.family=T" } else {  $booleanS .= " merge.family=F" }

		if ($parallel) {
			for my $chr (@chrl) {
				$chr="chr$chr";
				my $clusterf = "$cldir/$sample.$chr.cluster";
				my $job = "$sample.$rasym.$chr";
				my $logf = "$logdir/$sample.$rasym.log";
	    	if (! -e $clusterf) {
					print qq{
					echo "running $exec for $chr ..."
					bsub -R rusage[mem=10000] -q $q -W $time -g /$sample -J $job -o $logf \\
					"$R sample=\\\\\\"$sample\\\\\\\" dir=\\\\\\\"$bdir\\\\\\\" chr=\\\\\\\"$chr\\\\\\\" \\
					ref=\\\\\\\"$ref\\\\\\\" rasym=\\\\\\\"$rasym\\\\\\\" min.ram=$minram jittering=$jittering \\
					$booleanS < $tea/R/run.rid.r"
					};
	      } else {
		      print "echo \"skipping $exec for $chr: the cluster file already exists\"\n";
				}
			}
		} else {
	  		print qq{
			$R sample=\\"$sample\\" dir=\\"$bdir\\" ref=\\"$ref\\" rasym=\\"$rasym\\" \\
			min.ram=$minram jittering=$jittering $booleanS < $tea/R/run.rid.r
			};
		}
		#print qq {
		#	echo "pathogen enrichment..."
		#	python $tea/R/pathogen_enrich.py $bamdir/$sample $tea/lib/viruses.txt $dir/cluster_vam/$sample.enrich $dir/cluster_vam/$sample.cluster
		#}
	}
	if  ($step == 8)
	{
		print qq {
			echo "pathogen enrichment..."
			python $tea/R/pathogen_enrich.py $bamdir/$sample $tea/lib/viruses.txt $dir/cluster_vam/$sample.enrich $dir/cluster_vam/$sample.cluster
		}
	}

	# call germline
	if ($step <= 7){ 
	#if ($step <= 7 && $rasym ne "um") {
		my $R = "R --no-save --no-restore --no-environ --args ";
		if ($exec eq "rid.te.primate") {
			$R .= "te.primate.germline";
		} else {
			$R .= "germline";
		}
		my $booleanS;
		if ($oneside_ram) { $booleanS .= " oneside.ram=T" } else { $booleanS .= " oneside.ram=F" }
		if ($parallel) { $booleanS .= " parallel=T" } else { $booleanS .= " parallel=F" }
		if ($contig) { $booleanS .= " contig=T" } else { $booleanS .= " cointig=F" }
		my $min_out_conf = 5;
		##HC
		if ($exogenous) { $min_out_conf = 5 }; 

		if ($parallel) {
			my $job = "$sample.$rasym.germline";
			my $logf = "$logdir/$sample.$rasym.log";
			if ($step == 7) { 
			print qq{
				echo "calling germline events ..."
				bsub -q short -W 1:00 -g /$sample -J $job \\
				-o $logf \\
				"$R sample=\\\\\\\"$sample\\\\\\\" dir=\\\\\\\"$bdir\\\\\\\" ref=\\\\\\\"$ref\\\\\\\" \\
				rasym=\\\\\\\"$rasym\\\\\\\" min.ram=$minram min.out.conf=$min_out_conf $booleanS < $tea/R/run.rid.r"
			}
			} else {
			print qq{
				echo "calling germline events ..."
				bsub -q short -W 1:00 -g /$sample -J $job -w $sample.$rasym.chr* \\
				-o $logf \\
				"$R sample=\\\\\\\"$sample\\\\\\\" dir=\\\\\\\"$bdir\\\\\\\" ref=\\\\\\\"$ref\\\\\\\" \\
				rasym=\\\\\\\"$rasym\\\\\\\" min.ram=$minram min.out.conf=$min_out_conf $booleanS < $tea/R/run.rid.r"
			}
			}
		} else {
	  		print qq{
			$R sample=\\"$sample\\" dir=\\"$bdir\\" ref=\\"$ref\\" rasym=\\"$rasym\\" \\
			min.ram=$minram min.out.conf=$min_out_conf $booleanS < $tea/R/run.rid.r
			};
			if ($transduction)  {
					print qq{
						mv $bdir/$sample/cluster_umm $bdir/$sample/cluster_transduction
						python $tea/scripts/python/transduction_contig_module2.py $bdir/ $sample $bedtools $samtools $bwa $ra_um $reference 
					}
				} 
			if ($orphan)  {
					print qq{
						#mv $bdir/$sample/cluster_umm $bdir/$sample/cluster_orphan
						python $tea/scripts/python/orphan_contig_module2.py $bdir/ $sample $bedtools $samtools $bwa $ra_um $reference
						python $tea/scripts/python/post_process.py $bdir/ $sample $reference $tea $ra_um $bwa
					}
				} 
		}
	} 

	print "echo \"######## complete ##########\"\n";
}

sub bam_chr
{
    my $bamf = @_[0];
    my @headers = `$samtools view -H $bamf`;
    my $bam_chr = 0;
    for my $h (@headers) {
        if ($h =~ m/^\@SQ/) {
            $h =~ s/.*SN:(.*)\s.*/$1/;
            if ($h =~ m/chr/) {
                $bam_chr = 1;
                last;
            }
        }
    }
    return $bam_chr;
}

sub rmasker_coordinates
{
  my ($rmasker_file, $rmasker_start_line, $verbose) = @ARGV;
  my $outf = $rmasker_file;
  $outf =~ s/gz/coordinates.gz/;

  open(OUT, "| gzip -c > $outf") || die "can't create $outf: $!";

  # read rmasker and get the instance coordinates
  print "reading repeatmasker file $rmasker_file ...\n" if ($verbose);

  open(RM, "zcat $rmasker_file | ") || die "can't open repeat masker file $rmasker_file: $!";

  my $cnt=1;
  while ( $cnt++ < $rmasker_start_line) { <RM> }

  my $icnt = 0;
  my %reppos;
  my ($name, $divscore, $chr, $s, $e);

  while (<RM>) {
    chomp; my $line=$_;
    $line=~s/^\s+//;
    my @data = split(/\s+/,$line);
    #$name = "$data[9]:$data[10]:$data[11]";
		$data[10] =~ s/_$//;# remove the tailing _
    $name = "$data[10]:$data[11]:$data[12]";
    $divscore = $data[2]+$data[3]+$data[4];
    $chr = $data[5];
    $s = $data[6]+1; # note that all rmasker start coordinates are 0-t()based not 1-based
    $e = $data[7];

    push @{$reppos{uc($name)}}, "$chr:$s-$e";
    #push @{$reppos{$name}}, "$chr:$s-$e";
    $icnt++;
  }

 for my $n (keys %reppos) {
    print OUT "$n\t@{$reppos{$n}}\n";
  }
  close OUT;
  print "done extracting $icnt repeat instance coordinates of ".scalar(keys %reppos)." repeat types.\n" if ($verbose);
}

sub rmasker_coordinates_strand
{
  my ($rmasker_file, $rmasker_start_line, $verbose) = @ARGV;
  my $outf = $rmasker_file;
  $outf =~ s/gz/coordinates.strand.gz/;

  open(OUT, "| gzip -c > $outf") || die "can't create $outf: $!";

  # read rmasker and get the instance coordinates
  print "reading repeatmasker file $rmasker_file ...\n" if ($verbose);

  open(RM, "zcat $rmasker_file | ") || die "can't open repeat masker file $rmasker_file: $!";

  my $cnt=1;
  while ( $cnt++ < $rmasker_start_line) { <RM> }

  my $icnt = 0;
  my %reppos;
  my ($name, $divscore, $chr, $s, $e, $strand, $div);

  while (<RM>) {
    chomp; my $line=$_;
    $line=~s/^\s+//;
    my @data = split(/\s+/,$line);
    #$name = "$data[9]:$data[10]:$data[11]";
    $data[10] =~ s/_$//;# remove the tailing _
    $name = "$data[10]:$data[11]:$data[12]";
    $divscore = $data[2]+$data[3]+$data[4];
    $chr = $data[5];
    $s = $data[6]+1; # note that all rmasker start coordinates are 0-t()based not 1-based
    $e = $data[7];
	$strand = $data[9];

    push @{$reppos{uc($name)}}, "$chr:$s,$e:$strand:$divscore";
    #push @{$reppos{$name}}, "$chr:$s-$e";
    $icnt++;
  }

 for my $n (keys %reppos) {
    print OUT "$n\t@{$reppos{$n}}\n";
  }
  close OUT;
  print "done extracting $icnt repeat instance coordinates of ".scalar(keys %reppos)." repeat types.\n" if ($verbose);
}


sub make_subfamilyseq
{
  my ($family, $map_subfamily, $rdbs) = @_;
	my ($cmd, $db, $diag);
  if ($family eq "L1") { 
    $db = "/home/el114/repeat_analysis/data/consensus/M80343.1.fasta";
	} elsif (uc($family) eq "ALU") {
    $db = "/home/el114/repeat_analysis/data/consensus/AluY.fasta";
	} elsif ($family eq "PABL_A-INT") {
		$db = "/home/el114/repeat_analysis/data/consensus/PABL_A-INT.fasta"; 
	}
	push(@$rdbs, $db);
	if (!$map_subfamily) { return; }

  if ($family eq "L1") {
    open(L, "< $db") || die "can't open $db";
    my $seq = <L>; $seq = ""; while (<L>) { chomp; $seq = $seq.$_; } # get the consensus seq
    # print length($seq)."\n";

    $diag = "/home/el114/repeat_analysis/data/consensus/diagnostic_L1_positions.txt";
    open(D, "< $diag") || die "can't open $diag\n";
		my @m = <D>; chomp(@m);
		my @subfamily = split(/\t/, $m[0]); shift(@subfamily);
		my %seqs = (); foreach (@subfamily) { $seqs{$_} = $seq; } # subfamily consensus seq initialization
		my (@a, $pos);
		for (my $i=1; $i<@m; $i++) {
			@a = split(/\t/, $m[$i]);
			$pos = shift(@a);
			for (my $j=0; $j<@a; $j++) { if ($a[$j] ne "-") { 
				#print "replacing $subfamily[$j] $pos, to $a[$j]\n"; 
				substr($seqs{$subfamily[$j]}, $pos-1, 1, uc($a[$j])); } 
			}
		}
		foreach (@subfamily) {
				my $fasta = $db;
				$fasta =~ s/(.*).fasta/$1/;
				$fasta = "${fasta}_$_.fasta";	
				push(@$rdbs, $fasta);
				open(O, "> $fasta") || die "can't create $fasta";
				print O ">$_\n$seqs{$_}\n";
				$cmd = "$blast_path/makeblastdb -dbtype nucl -in $fasta";	
				print "generating blast db for $fasta..\n";
				system($cmd);
				print "done generating.\n";
		}
	} elsif (uc($family) eq "Alu") {
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluYa5.fasta"); 
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluYb8.fasta"); 
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluYb9.fasta"); 
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluSc.fasta"); 
	} 
}

# blast clipped and ram mate contigs to consensus sequences
# and annotate insertion size, inversion status and diagnostic features for L1
# requires unique ids, otherwise generates ids [family_##]
sub map_consensus
{
	my ($file, $family, $map_subfamily, $verbose) = @_;

	open(F, "< $file") || die "can'e open $file";
	open(O, "> $file.fasta") || die "can't create $file.fasta"; # fasta
	open(O2, "> $file.annot") || die "can't create $file.annot"; # full annot 
	my $header = <F>; chomp($header);
	my @a = split(/\t/, $header); my %idx = ();
	for (my $i=0; $i<@a; $i++) { $idx{$a[$i]} = $i; }
		
	my (%insertions, $sname, $seq); # id to record
	print "reading $file and generating $file.fasta..\n";
	my @f = <F>;
	my @ids = ();
	if (exists $idx{id}) {
	for (my $i=0; $i<@f; $i++) {
		chomp($f[$i]); my @a=split(/\t/, $f[$i]);
		if ($a[$idx{family}] eq $family) {
		push(@ids, $a[$idx{id}]);
		$insertions{$a[$idx{id}]}{line} = $f[$i];
		$insertions{$a[$idx{id}]}{orientation} = $a[$idx{orientation}]; 
		$insertions{$a[$idx{id}]}{polyA} = $a[$idx{polyA}]; 
		$insertions{$a[$idx{id}]}{polyT} = $a[$idx{polyT}]; 

		$sname = "$a[$idx{id}]:pclipped";
		$seq = $a[$idx{pclipped}];
		print O ">$sname\n$seq\n";

		$sname = "$a[$idx{id}]:nclipped";
		$seq = $a[$idx{nclipped}];
		print O ">$sname\n$seq\n";

		$sname = "$a[$idx{id}]:prammate";
		$seq = $a[$idx{prammate}];
		print O ">$sname\n$seq\n";

		$sname = "$a[$idx{id}]:nrammate";
		$seq = $a[$idx{nrammate}];
		print O ">$sname\n$seq\n";
	}
	}
	} else { # generate unique id
	my ($cnt, $uid) = (1, "");
	for (my $i=0; $i<@f; $i++) {
		chomp($f[$i]); my @a=split(/\t/, $f[$i]);
		# for the modification for the young repeat families (L1HS) and polyA support
		#if ($a[$idx{family}] eq $family) {
		if ($a[$idx{family}] =~ m/$family/ || $a[$idx{"rep.repeat"}] =~ m/$family/) {
		$uid = "$family\_$cnt"; $cnt++;
		push(@ids, $uid);
		$insertions{$uid}{line} = "$uid\t$f[$i]";
		$insertions{$uid}{orientation} = $a[$idx{orientation}]; 
		$insertions{$uid}{polyA} = $a[$idx{polyA}]; 
		$insertions{$uid}{polyT} = $a[$idx{polyT}]; 

		$sname = "$uid:pclipped";
		$seq = $a[$idx{pclipped}];
		print O ">$sname\n$seq\n";

		$sname = "$uid:nclipped";
		$seq = $a[$idx{nclipped}];
		print O ">$sname\n$seq\n";

		$sname = "$uid:prammate";
		$seq = $a[$idx{prammate}];
		print O ">$sname\n$seq\n";
	}
	}
		$header = "id\t$header";
	}
	print "done generating.\n";
	close F; close O;
	
	my @dbs = (); # db for each subfamily
	&make_subfamilyseq($family, $map_subfamily, \@dbs);
	my @asf = ();
	for my $sf (@dbs) {
		my $subfamily = basename($sf); $subfamily =~ s/(.*).fasta/$1/;
		push(@asf, $subfamily);
		print "subfamily: $subfamily\n";
		print "mapping consensus for $subfamily\n";
		&map_consensus_db(\%insertions, "$file.fasta", $sf, $subfamily, $verbose); # inversion, insertion size, %divergence 
	}
	print O2 "$header\tcontig_orientation\torientation\tinversion\tsize\tcstart\tcend\tsubfamily\tmismatch\tdivpct\n";
	for my $id (@ids) {
		print O2 "$insertions{$id}{line}";

		# summarize divs and extract subfamilies with the minimum div		
		my (@divs, @divpcts, @est_subfamily);
		my $mindiv = 9999;
		for my $sf (@asf) {
			my $ssf = $sf; # short sf removing M80343.1 prefix
			#if (($family eq "L1" || $family eq "L1HS" ) && $sf ne "M80343.1") { # skip L1 ref name
			if ($family =~  m/L1/ && $sf ne "M80343.1") { # skip L1 ref name
				$ssf =~ s/M80343.1_//;
			}
			push(@divs, "$ssf:$insertions{$id}{$sf}{div}");
			push(@divpcts, "$ssf:$insertions{$id}{$sf}{divpct}");
			if ($insertions{$id}{$sf}{div} < $mindiv) { $mindiv = $insertions{$id}{$sf}{div} }
		}
		for my $sf (@asf) {
			if ($sf eq "M80343.1") { next } 
			my $ssf = $sf; $ssf =~ s/M80343.1_//;
			if ($insertions{$id}{$sf}{div} == $mindiv) { push(@est_subfamily, $ssf) }
		}
		if (scalar @est_subfamily == (scalar @divs)-1) { @est_subfamily = (); }  
	
		# decide the final orientation from polyA/T and contig mapping	
		my $orientation = $insertions{$id}{orientation};
		if ($orientation eq "NA" && exists $insertions{$id}{$asf[0]}{map_orientation}) { 
			$orientation = $insertions{$id}{$asf[0]}{map_orientation};
		}
		if (exists $insertions{$id}{$asf[0]}{map_orientation}) {
			print O2 "\t$insertions{$id}{$asf[0]}{map_orientation}\t$orientation\t$insertions{$id}{$asf[0]}{inv}\t$insertions{$id}{$asf[0]}{size}\t$insertions{$id}{$asf[0]}{cstart}\t$insertions{$id}{$asf[0]}{cend}\t";
			if (@est_subfamily >=1) {
				print O2 join(",", @est_subfamily)."\t";
			} else {
				print O2 "-\t"; 
			}
			print O2 join(",", @divs)."\t";
			print O2 join(",", @divpcts)."\n";
		} else {
			print O2 "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
		}
	}
	close O2;
	`rm $file.fasta`;
	if (!$verbose) { `rm $file.*.fasta.map` }
}

sub map_consensus_db
{
	my ($href, $fasta, $db, $subfamily, $verbose) = @_;

	print "blasting $fasta to $db..\n";
	my $cmd = "$blast_path/legacy_blast.pl blastall -p blastn -m 8 -i $fasta -d $db --path $blast_path > $fasta.".basename($db).".map";
	print $cmd."\n"; system($cmd);
	print "done.\n";

	open(F, "< $db") || die "can't open $db";
	my $l = <F>; $l = ""; while (<F>) { chomp; $l = $l.$_; } $l = length($l);
	&parse_map($href, $fasta.".".basename($db).".map", $l, $subfamily, $verbose);	
	if (!$verbose) {
		$cmd = "rm $fasta.".basename($db).".map";
		#system($cmd);
	}
}

sub parse_map
{
	my ($href, $mapfile, $full_length, $sf, $verbose) = @_;
	# fille the hash for each id; if two mappings then chose the longer mapping
	open(F, "< $mapfile") || die "can't open $mapfile";
	print "reading $mapfile\n";

	while (<F>) {
		chomp;
		my ($qid, $sid, $identity, $alnlen, $mismatch, $gap, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/);
		#print "$qid, $sid, $identity, $alnlen, $mismatch, $gap, $qstart, $qend, $sstart, $send, $evalue, $bitscore\n";
		my ($id, $type) = ($qid, $qid); $id =~ s/(.*):.*/$1/; $type =~ s/.*:(.*)/$1/; 
		if (!exists$href->{$id}{$sf}{$type} || ($alnlen > $href->{$id}{$sf}{$type}{alnlen}) || ($$alnlen == $href->{$id}{$sf}{$type}{alnlen} && $href->{$id}{$sf}{$type}{mismatch} > $mismatch))  { # choose the maximum match and minimum mismaches where thee are multiple mappings
			$href->{$id}{$sf}{$type}{alnlen} = $alnlen; 
			$href->{$id}{$sf}{$type}{mismatch} = $mismatch; 
			$href->{$id}{$sf}{$type}{gap} = $gap; 
			$href->{$id}{$sf}{$type}{sstart} = $sstart; 
			$href->{$id}{$sf}{$type}{send} = $send; 
			if ($sstart <= $send) {
				$href->{$id}{$sf}{$type}{orientation} = "+"; 
			} else { 
				$href->{$id}{$sf}{$type}{orientation} = "-"; 
			}
		}
	}

	for my $id (keys %{$href}) {
		my ($orientation, $inv, $start, $end, $size, $alnlen, $div, $divpct) = ("NA", "NA", -1, -1, -1, 0, 0, -1);
		if (exists $href->{$id}{$sf}) {
		for my $type (keys %{$href->{$id}{$sf}}) {
			for my $key (keys %{$href->{$id}{$sf}{$type}}) {
				if ($verbose) { print "$id->$sf->$type->$key: $href->{$id}{$sf}{$type}{$key}\n"; }
			}
			$div = $div + $href->{$id}{$sf}{$type}{mismatch} + $href->{$id}{$sf}{$type}{gap};
			$alnlen =  $alnlen + $href->{$id}{$sf}{$type}{alnlen};
		}
		# orientation: + or - if both prammate (nclipped)  and nrammate (pclipped) are + or -
		my ($l, $r) = ("NA", "NA");
		if (exists $href->{$id}{$sf}{prammate}{orientation}) { 
			$l = $href->{$id}{$sf}{prammate}{orientation}; 
		} elsif (exists $href->{$id}{$sf}{nclipped}{orientation}) { 
      $l = $href->{$id}{$sf}{nclipped}{orientation};
		} 
		if (exists $href->{$id}{$sf}{nrammate}{orientation}) { 
			$r = $href->{$id}{$sf}{nrammate}{orientation}; 
		} elsif (exists $href->{$id}{$sf}{pclipped}{orientation}) { 
      $r = $href->{$id}{$sf}{pclipped}{orientation};
		}
		if ($l ne "NA" && $r ne "NA") {
			if ($l eq "+" && $r eq "+") {
				$orientation = "+"; 
				$start = 9999;	
				if (exists $href->{$id}{$sf}{prammate}{sstart} && $href->{$id}{$sf}{prammate}{sstart}< $start) { $start = $href->{$id}{$sf}{prammate}{sstart}}
				if (exists $href->{$id}{$sf}{nclipped}{sstart} && $href->{$id}{$sf}{nclipped}{sstart}< $start) { $start = $href->{$id}{$sf}{nclipped}{sstart}}
				if (exists $href->{$id}{$sf}{nrammate}{sstart} && $href->{$id}{$sf}{nrammate}{sstart}< $start) { $start = $href->{$id}{$sf}{nrammate}{sstart}}
				$end = -1;
				if (exists $href->{$id}{$sf}{nrammate}{send} && $href->{$id}{$sf}{nrammate}{send} > $end) { $end = $href->{$id}{$sf}{nrammate}{send}}
				if (exists $href->{$id}{$sf}{pclipped}{send} && $href->{$id}{$sf}{pclipped}{send} > $end) { $end = $href->{$id}{$sf}{pclipped}{send}}
				if (exists $href->{$id}{$sf}{prammate}{send} && $href->{$id}{$sf}{prammate}{send} > $end) { $end = $href->{$id}{$sf}{prammate}{send}}
				$size = $end - $start + 1;
				$inv = "NO_INV";
			} elsif ($l eq "-" && $r eq "-") {
				$orientation = "-"; 
				$start = 9999;
				if (exists $href->{$id}{$sf}{nrammate}{send} && $href->{$id}{$sf}{nrammate}{send}< $start) { $start = $href->{$id}{$sf}{nrammate}{send}}
				if (exists $href->{$id}{$sf}{pclipped}{send} && $href->{$id}{$sf}{pclipped}{send}< $start) { $start = $href->{$id}{$sf}{pclipped}{send}}
				if (exists $href->{$id}{$sf}{prammate}{send} && $href->{$id}{$sf}{prammate}{send}< $start) { $start = $href->{$id}{$sf}{prammate}{send}}
				$end = -1;
				if (exists $href->{$id}{$sf}{prammate}{sstart} && $href->{$id}{$sf}{prammate}{sstart} > $end) { $end = $href->{$id}{$sf}{prammate}{sstart}}
				if (exists $href->{$id}{$sf}{nclipped}{sstart} && $href->{$id}{$sf}{nclipped}{sstart} > $end) { $end = $href->{$id}{$sf}{nclipped}{sstart}}
				if (exists $href->{$id}{$sf}{nrammate}{sstart} && $href->{$id}{$sf}{nrammate}{sstart} > $end) { $end = $href->{$id}{$sf}{nrammate}{sstart}}
				$inv = "NO_INV";
			} else {
				$inv = "INV";
			}
		} elsif ($l eq "NA" && $href->{$id}{polyT} eq "polyT") {
			$orientation = "-";
			if ($r eq "-") {
				$inv = "NO_INV";
				$start = 9999; 
				if (exists $href->{$id}{$sf}{nrammate}{send} && $href->{$id}{$sf}{nrammate}{send}< $start) { $start = $href->{$id}{$sf}{nrammate}{send}}
				if (exists $href->{$id}{$sf}{pclipped}{send} && $href->{$id}{$sf}{pclipped}{send}< $start) { $start = $href->{$id}{$sf}{pclipped}{send}}
				if (exists $href->{$id}{$sf}{prammate}{send} && $href->{$id}{$sf}{prammate}{send}< $start) { $start = $href->{$id}{$sf}{prammate}{send}}
				$end = $full_length;
			} else {
				$inv = "INV";
			}
		} elsif ($r eq "NA" && $href->{$id}{polyA} eq "polyA") {
			$orientation = "+";
			if ($l eq "+") {
				$inv = "NO_INV";
				$start = 9999;
				if (exists $href->{$id}{$sf}{prammate}{sstart} && $href->{$id}{$sf}{prammate}{sstart}< $start) { $start = $href->{$id}{$sf}{prammate}{sstart}}
				if (exists $href->{$id}{$sf}{nclipped}{sstart} && $href->{$id}{$sf}{nclipped}{sstart}< $start) { $start = $href->{$id}{$sf}{nclipped}{sstart}}
				if (exists $href->{$id}{$sf}{nrammate}{sstart} && $href->{$id}{$sf}{nrammate}{sstart}< $start) { $start = $href->{$id}{$sf}{nrammate}{sstart}}
				$end = $full_length;
			} else {
				$inv = "INV";
			}
		} 
			if ($start != -1 && $end != -1) {	$size = $end - $start + 1; }
			$divpct = sprintf("%.1f", $div/$alnlen*100); 
			$href->{$id}{$sf}{map_orientation} = $orientation;
			$href->{$id}{$sf}{inv} = $inv;
			$href->{$id}{$sf}{cstart} = $start;
			$href->{$id}{$sf}{cend} = $end;
			$href->{$id}{$sf}{size} = $size;
			$href->{$id}{$sf}{divpct} = $divpct;
			$href->{$id}{$sf}{div} =  $div; 
			print "$l: $r: $id: $sf: $orientation, $inv, $size, $start, $end, $div, $divpct, $alnlen\n";
	}
	}
}

sub contig {

	my $usage = qq {
	Usage:   contig [options] <input file>
 
	Options: -a STR path to the contig assembler (cap3) 
	         -d STR path to RUN DIRECTORY (upper to the sample folder) 
	              with read names and clipped sequences and disc. bam files)
	         -o STR outdir (.)
	         -C INT confidence level for filtering (0)
	         -m STR consensus mapping [L1|Alu|SVA] (default: no)
	         -s subfamily mapping (default: no)
	         -t STR tempdir for running cap3 (/hms/scratch1/el114)
	         -c     clipped sequence contig only (0)
	         -r     ram mate contig only (0)
	         -e     extended ram mate contig only (0)
	         -x     exogenous (0)		
	         -p     Tea was run for each chromosome (0)		
	         -R STR  assembly symbol
	         -v     verbose

	};

	my %opts = ();
	getopts("ha:d:m:o:vrcevst:C:xR:p", \%opts);

	if (@ARGV < 1 || exists($opts{h})) { die $usage };
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
	my $conly = 0; $conly = 1 if (defined $opts{c});
	my $ronly = 0; $ronly = 1 if (defined $opts{r});
	my $eonly = 0; $eonly = 1 if (defined $opts{e});
	my $map = "no"; $map = $opts{m} if (defined $opts{m});
	my $map_subfamily = 0; $map_subfamily = 1 if (defined $opts{s});
	my $assembler = $assembler; $assembler = $opts{a} if (defined $opts{a});
	my $dir = ""; $dir = $opts{d} if (defined $opts{d});
	my $odir = "."; $odir = $opts{o} if (defined $opts{o});
	my $tdir = "$temp_loc/$map"; $tdir = $opts{t} if (defined $opts{t});
	my $conf = 0;  $conf = $opts{C} if (defined $opts{C}); 
	if (! -d $tdir) { system("mkdir -p $tdir") }
	my $exogenous = 0; $exogenous = 1 if (defined $opts{x});
	my $rasym = ""; $rasym = $opts{R} if (defined($opts{R}));
	my $parallel = 0; $parallel = 1 if (defined($opts{p}));

	my $infile = shift(@ARGV);

	my $outfile = $odir."/".basename($infile).".contig";

	open(F, "< $infile") || die "can't open $infile";

	my %can = (); # a hash for the candidate regions
	my @can = <F>; chomp(@can);
	close F;

	my $header = $can[0]; # record the column index
	my @a = split(/\t/, $header); my %idx = ();
	for (my $i=0; $i<@a; $i++) {
		$idx{$a[$i]} = $i;
	}

	# if there is no sample column, extract from the input filename by dropping .cluster or .germline or .somatic
	my $sample = "";
	if (!exists $idx{sample}) {
		$sample = basename($infile);	
		$sample =~ s/(.cluster|.germline|.somatic).*//;
	}

	for (my $i=1; $i<@can; $i++) {
		my @a = split(/\t/, $can[$i]);
		if ($sample eq "") { $sample = $a[$idx{sample}] }; # if sample column exists
		if ($conf>0) { # filter candidates with the confidence level lower than the specified
			if (!exists $idx{conf}) { 
				die("no conf column in the input file")
			} elsif ($a[$idx{conf}] <$conf) { next }
		}
		if (!exists $can{$sample}) {$can{$sample}{chr} = () }
		if (!exists $can{$sample}{$a[$idx{chr}]}) { 
			push(@{$can{$sample}{chr}}, $a[$idx{chr}]);
			$can{$sample}{$a[$idx{chr}]}{site} = ();
		}

		my ($site, $pramsite, $nramsite);
		if (exists $idx{s} && exists $idx{e}) {
			$site = "$a[$idx{s}]:$a[$idx{e}]";
		} else {
			# for the compatibility with the old column format
			$site = "$a[$idx{pram_start}]:$a[$idx{nram_end}]";
		}
		$pramsite = "$a[$idx{pram_start}]-$a[$idx{pram_end}]";
		$nramsite = "$a[$idx{nram_start}]-$a[$idx{nram_end}]";
  		push(@{$can{$sample}{$a[$idx{chr}]}{site}}, $site);
  		$can{$sample}{$a[$idx{chr}]}{$site}{pramsite} = $pramsite;
  		$can{$sample}{$a[$idx{chr}]}{$site}{nramsite} = $nramsite;
  		$can{$sample}{$a[$idx{chr}]}{$site}{type} = $a[$idx{"rep.repeat"}];
  		$can{$sample}{$a[$idx{chr}]}{$site}{rec} = $can[$i];
	}

	for my $s (keys %can) {
		print "processing $s\n" if ($verbose);
	    &clipped_contig($s, \%can, $dir, $rasym, $parallel, $verbose) if (!$ronly && !$eonly);
	    &ram_contig($s, \%can, $dir, $exogenous, $rasym, $parallel, $verbose) if (!$conly && !$eonly);
	    &eram_contig($s, \%can, $dir, $exogenous, $verbose) if (!$conly && !$ronly);
	} 

	print "generating the contigs...\n" if ($verbose);
	&make_contig2(\%can, $odir, $tdir, $ronly, $conly, $eonly, $assembler, $verbose);
	print "... done\n" if ($verbose);

	print "writing the contigs...\n" if ($verbose);
	&write_contig2($outfile, \%can, $header);
	print "... done\n" if ($verbose);

	if ($map ne "no") { 
		print "mapping the contigs to the consensus sequences... \n" if ($verbose);
		&map_consensus($outfile, $map, $map_subfamily, $verbose); 
		print "...done\n" if ($verbose);
	} 
}

sub ram_contig_line 
{
	my ($s, $href, $fh, $chr, $hprn, $hnrn, $hpcrn, $hncrn, $verbose) = @_;
	while (<$fh>) {
		chomp; my @a = split(/\t/);
		$chr = $a[0];
		if ($chr eq "") { 
			$chr = $a[0]; 
			#unless ($chr =~ m/^chr/) { $chr = "chr$chr"; }; 
		}
		#$chr = "chr$chr";
        my $hhref = $href->{$s}{$chr};
        if (exists $hhref->{"$a[1]:$a[2]"}) {
            my @pr = split(",", $a[20]);
            my @nr = split(",", $a[21]);
            my @ppos = split(",", $a[14]);
            my @npos = split(",", $a[15]);
            my $region = "$chr:$a[1]:$a[2]";
            for (my $i=0; $i<@pr; $i++) {
                if ($pr[$i] =~ m/mu1$/ || $pr[$i] =~ m/mu2$/ || $pr[$i] =~ m/sc$/) {
                    $hpcrn->{$pr[$i]}{chr} = $chr;
                    $hpcrn->{$pr[$i]}{pos} = $ppos[$i]; # we need to look the seq and qual for the reads mapped to the neg
                    # Note that one read can be mapped to multiple clusters!
                    push(@{$hpcrn->{$pr[$i]}{region}}, $region);
                } else {
                    $hprn->{$pr[$i]}{chr} = $chr;
                    $hprn->{$pr[$i]}{pos} = $ppos[$i];
                    push(@{$hprn->{$pr[$i]}{region}}, $region);
                }    
            }    
            for (my $i=0; $i<@nr; $i++) {
                if ($nr[$i] =~ m/mu1$/ || $nr[$i] =~ m/mu2$/ || $nr[$i] =~ m/sc$/) {
                    $hncrn->{$nr[$i]}{chr} = $chr;
                    $hncrn->{$nr[$i]}{pos} = $npos[$i];
                    push(@{$hncrn->{$nr[$i]}{region}}, $region);
                } else {
                    $hnrn->{$nr[$i]}{chr} = $chr;
                    $hnrn->{$nr[$i]}{pos} = $npos[$i];
                    push(@{$hnrn->{$nr[$i]}{region}}, $region);
                }    
            }    
            push(@{$hhref->{"$a[1]:$a[2]"}->{pr}}, @pr); # positive ram repeat reads
            push(@{$hhref->{"$a[1]:$a[2]"}->{nr}}, @nr); # netagive ram repeat reads
            push(@{$hhref->{"$a[1]:$a[2]"}->{ppos}}, @ppos); # positive ram repeat reads
            push(@{$hhref->{"$a[1]:$a[2]"}->{npos}}, @npos); # netagive ram repeat reads

			if ($verbose) {
				print "$s:$chr:$a[1]:$a[2]:".join(",",@{$hhref->{"$a[1]:$a[2]"}->{pr}})."\n";
				print "$s:$chr:$a[1]:$a[2]:".join(",", @{$hhref->{"$a[1]:$a[2]"}->{nr}})."\n";
				print "$s:$chr:$a[1]:$a[2]:".join(",",@{$hhref->{"$a[1]:$a[2]"}->{ppos}})."\n";
				print "$s:$chr:$a[1]:$a[2]:".join(",", @{$hhref->{"$a[1]:$a[2]"}->{npos}})."\n";
			}
        }
    }
}

sub eram_contig_line 
{
	my ($s, $href, $fh, $chr, $site, $orientation, $hr, $verbose) = @_;
	$site =~ s/-/:/; 
	my $region = $site; $region = "$chr:$region";
    my $hhref = $href->{$s}{$chr};
	while (<$fh>) {
		chomp; my @a = split(/\t/);
		unless ($chr =~ m/^chr/) { $chr = "chr$chr"; }; 
		if ($orientation == 0) { # nram
			push(@{$hhref->{$site}->{enr}}, $a[0]); 
			push(@{$hhref->{$site}->{enpos}}, $a[3]);
		} else { # pram
			push(@{$hhref->{$site}->{epr}}, $a[0]); 
			push(@{$hhref->{$site}->{eppos}}, $a[3]); 
		}

		$hr->{$a[0]}{chr} = $chr;
		$hr->{$a[0]}{pos} = $a[3];
		push(@{$hr->{$a[0]}{region}}, $region);
	}
}

sub eram_contig
{
    my ($s, $href, $dir, $exogenous, $verbose) = @_;

    my $bamdir = "$dir/$s/bam";
	my ($bamf, $clbamf, $umbamf) = ("$bamdir/$s.disc.bam", "$bamdir/$s.cl.disc.bam", "$bamdir/$s.um.bam");
    my (%prn, %nrn, %pcrn, %ncrn); # caution: one read can be mapped to multiple clusters!
	my $bam_chr = 1;
	if ($exogenous) { $bam_chr = &bam_chr($umbamf); } else { $bam_chr = &bam_chr($bamf); };

	for my $chr (@{$href->{$s}{chr}}) { # for each chromosome
		my $bchr = $chr; if ($bam_chr == 0) { $bchr =~ s/chr//; }
		for my $site (@{$href->{$s}{$chr}{site}}) { # for each candidate site
			my ($pramsite, $nramsite) = ($href->{$s}{$chr}{$site}{pramsite}, $href->{$s}{$chr}{$site}{nramsite});
			if ($verbose) { print "site: $chr:$site:$pramsite:$nramsite\n"; }
			if ($exogenous == 0) { 
				# check if the below checking is enough
				#$bchr =~ s/chr//;
				if ($pramsite ne "0-0") {
				
				open(my $fh, "$samtools view -F 0x4 -F 0x10 $bamf $bchr:$pramsite |");
				&eram_contig_line($s, $href, $fh, $chr, $site, 1, \%prn, $verbose);
				close $fh;

				open(my $fh, "$samtools view -F 0x4 -F 0x10 $clbamf $bchr:$pramsite |");	
				&eram_contig_line($s, $href, $fh, $chr, $site, 1, \%pcrn, $verbose);
				close $fh;
				}

				if ($nramsite ne "0-0") {
                open(my $fh, "$samtools view -F 0x4 -f 0x10 $bamf $bchr:$nramsite |");
                &eram_contig_line($s, $href, $fh, $chr, $site, 0, \%nrn, $verbose);
                close $fh;

                open(my $fh, "$samtools view -F 0x4 -f 0x10 $clbamf $bchr:$nramsite |");
                &eram_contig_line($s, $href, $fh, $chr, $site, 0, \%ncrn, $verbose);
                close $fh;
				}
			} else {
				if ($pramsite ne "0-0") {
				open(my $fh, "$samtools view -F 0x4 -F 0x10 $umbamf $bchr:$pramsite |");	
				&eram_contig_line($s, $href, $fh, $chr, $site, 1, \%pcrn, $verbose);
				close $fh;
				}

				if ($nramsite ne "0-0") {
                open(my $fh, "$samtools view -F 0x4 -f 0x10 $umbamf $bchr:$nramsite |");
                &eram_contig_line($s, $href, $fh, $chr, $site, 0, \%ncrn, $verbose);
                close $fh;
				}
			}
		}
    }

	# read a disc or cl.disc.bam to get the read seq and qual
	if ($exogenous == 0) {
		print "looking for ".scalar keys(%prn).":".scalar keys(%nrn)." ram names from disc.bam\n";
		print "looking for ".scalar keys(%pcrn).":".scalar keys(%ncrn)." ram names from cl.disc.bam\n";

		&get_read_seq($href->{$s}, $bamf, \%prn, \%nrn, 1, $verbose);
		&get_read_seq($href->{$s}, $clbamf, \%pcrn, \%ncrn, 1, $verbose);
	} else {
		print "looking for ".scalar keys(%pcrn).":".scalar keys(%ncrn)." ram names from um.bam\n";
		&get_read_seq($href->{$s}, $umbamf, \%pcrn, \%ncrn, 1, $verbose);
	}
}

sub ram_contig
{
    my ($s, $href, $dir, $exogenous, $rasym, $parallel, $verbose) = @_;

    my $cldir;
    if ($rasym eq "") { $cldir = "$dir/$s/cluster" } else { $cldir = "$dir/$s/cluster_${rasym}m" }

    # extract ram read names for candidate regions
    my (%prn, %nrn, %pcrn, %ncrn); # caution: one read can be mapped to multiple clusters!

    if (!$parallel) {
        my $file = "$cldir/$s.cluster.raw";
        if (! -e $file) { print "$s does not have the cluster raw file $file!\n"; next; }
        open(my $fh, "< $file") || die "can't open $file";
        print "opened $file\n";
        <$fh>;  # discard the header
		&ram_contig_line($s, $href, $fh, "", \%prn, \%nrn, \%pcrn, \%ncrn, $verbose);
    } else {
    	for my $chr (@{$href->{$s}{chr}}) {
       		# get the ram read names
        	my $file = "$cldir/$s.$chr.cluster.raw";
        	if (! -e $file) { print "$s $chr does not have the cluster raw file!\n"; next; }
        	open(my $fh, "< $file") || die "can't open $file";
        	print "opened $file\n";
			<$fh>
			&ram_contig_line($s, $href, $fh, $chr, \%prn, \%nrn, \%pcrn, \%ncrn, $verbose);
    	}
    }
   # read a disc or cl.disc.bam to get the read seq and qual
    if ($exogenous == 0) {
        my ($bamf, $clbamf);
        $bamf = "$dir/$s/bam/$s.disc.bam";
        $clbamf = "$dir/$s/bam/$s.cl.disc.bam";
        print "looking for ".scalar keys(%prn).":".scalar keys(%nrn)." ram names from disc.bam\n";
        print "looking for ".scalar keys(%pcrn).":".scalar keys(%ncrn)." ram names from cl.disc.bam\n";

        &get_read_seq($href->{$s}, $bamf, \%prn, \%nrn, 0, $verbose);
        &get_read_seq($href->{$s}, $clbamf, \%pcrn, \%ncrn, 0, $verbose);
    } else {
        my $umbamf = "$dir/$s/bam/$s.um.bam";
        print "looking for ".scalar keys(%pcrn).":".scalar keys(%ncrn)." ram names from um.bam\n";

        &get_read_seq($href->{$s}, $umbamf, \%pcrn, \%ncrn, 0, $verbose);
    }
}


sub get_read_seq
{
	my ($href, $bamf, $phr, $nhr, $ext, $verbose) = @_; # if $ext == 1, extended ram contig
	my $x = scalar keys %{$phr};
	my $y = scalar keys %{$nhr};
	if ($x == 0 && $y == 0) { 
		print "no read to search\n";	
		return(); 
	}	
	my ($seq, $qual);
	my $um = 0; if ($bamf =~ /um.bam/) { $um = 1 }; # search an unmapped bam
	print "start reading $bamf for assembling repeat read sequences\n";
	open(F, "$samtools view -X $bamf |");
	my (@a, $n, $s, @b);
	while (<F>) { 
		if ($x == 0 && $y == 0) { last; }
		@a = split(/\t/); $n = $a[0]; # read name
		if (exists $phr->{$a[0]}) { # if a read name is in the ram list
			my $chr = $a[2];
			if (!($chr =~ m/chr/)) { $chr = "chr$chr" }; # for the bam without chr prefix
			if ($verbose) { print "$phr->{$a[0]}->{chr} <=> $chr, ".abs($phr->{$a[0]}->{pos})." <=> $a[3]\n"; }
			if (($um == 1 && $a[1] =~ /u/) || ($um == 0 && ($phr->{$a[0]}->{chr} ne $chr || abs($phr->{$a[0]}->{pos}) != $a[3]))) {
			# pos ram mate
			for my $r (@{$phr->{$a[0]}->{region}})	{
				my @b = split(":", $r);
				#if ($um == 0 && ! ($a[1] =~ m/r/))  {
				if (! ($a[1] =~ m/r/))  { # Sep.11, 2013 I'm not sure why I had $um=0 as above ; note that I always write the mate of pram in negative orientation
                                          # to make the notation consistent with nclipped
    				($a[9] = reverse $a[9]) =~ tr/gatcGATC/ctagCTAG/; # if flag has 'r' reverse complement
    				$a[10] = reverse $a[10];
				}
				if ($verbose) { print "prseq: found seq for $r, $a[0], $a[9]\n"; }
				if ($ext) {
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{eprseq}}, $a[9]);
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{eprqual}}, $a[10]);
				} else {
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{prseq}}, $a[9]);
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{prqual}}, $a[10]);
				}
			}
			delete($phr->{$a[0]});
			$x = scalar keys %{$phr};
			}
		}
		if (exists $nhr->{$a[0]}) {
			my $chr = $a[2];
			if (!($chr =~ m/chr/)) { $chr = "chr$chr" }; # for the bam without chr prefix
			if (($um == 1 && $a[1] =~ /u/) || ($um == 0 && ($nhr->{$a[0]}->{chr} ne $chr || abs($nhr->{$a[0]}->{pos}) != $a[3]))) {
			for my $r (@{$nhr->{$a[0]}->{region}})	{
				my @b = split(":", $r);
				#if ($um == 0 && $a[1] =~ m/r/)  {
				if ($a[1] =~ m/r/)  {
    			($a[9] = reverse $a[9]) =~ tr/gatcGATC/ctagCTAG/;
    			$a[10] = reverse $a[10];
				}
				if ($verbose) {print "nrseq: found seq for $r, $a[0], $a[9]\n"; }
				if ($ext) {
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{enrseq}}, $a[9]);
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{enrqual}}, $a[10]);
				} else {
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{nrseq}}, $a[9]);
					push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{nrqual}}, $a[10]);
				}
			}
			delete($nhr->{$a[0]});
			$y = scalar keys %{$nhr};
			}
		}
	}
	print "done reading $bamf for assembling repeat read sequences\n";
	close F;
}

# load clipped sequences and qual scores of clipped reads into a hash
sub clipped_contig
{
	my ($s, $href, $dir, $rasym, $parallel, $verbose) = @_;

	my $cldir;
	if ($rasym eq "") { $cldir = "$dir/$s/cluster" } else { $cldir = "$dir/$s/cluster_${rasym}m" }

	if (!$parallel) {
		my $clfile = "$cldir/$s.clipped";
		open(F, "< $clfile") || die "can't open $clfile";
		print "opened $clfile\n";
		<F>;
		while (<F>) {
			chomp; my @a = split(/\t/);
			my $chr = $a[0]; unless ($chr =~ m/^chr/) { $chr = "chr$chr"; };
			if (exists $href->{$s}->{$chr}->{"$a[1]:$a[2]"}) {
				if ($verbose) {	print "$s:$chr:$a[1]:$a[2]:"; }
				if ($a[3]>0) {
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pseq}}, $a[8]);
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pqual}}, $a[9]);
					if ($verbose) {	print "pseq".join(",", @{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pseq}})."\n"; }
				} else {
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nseq}}, $a[8]);
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nqual}}, $a[9]);
					if ($verbose) {	print "nseq".join(",", @{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nseq}})."\n"; }
				}
			}
		}
	} else {
	# get clipped sequences and read names for each candidate region
	for my $chr (@{$href->{$s}{chr}}) {
		my $clfile = "$cldir/$s.$chr.clipped";

		if (! -e $clfile) { print "$s $chr does not have the clipped sequence file $clfile!\n"; next; }
		open(F, "< $clfile") || die "can't open $clfile";
		print "opened $clfile\n";
		while (<F>) {
			chomp; my @a = split(/\t/);
			if (exists $href->{$s}->{$chr}->{"$a[1]:$a[2]"}) {
				if ($verbose) {	print "$s:$chr:$a[1]:$a[2]:"; }
				if ($a[3]>0) {
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pseq}}, $a[8]);
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pqual}}, $a[9]);
					if ($verbose) {	print "pseq".join(",", @{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pseq}})."\n"; }
				} else {
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nseq}}, $a[8]);
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nqual}}, $a[9]);
					if ($verbose) {	print "nseq".join(",", @{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nseq}})."\n"; }
				}
			}
		}
	}
	}
}

sub write_contig2
{
  my ($outfile, $href, $header) = @_;
 
  open(F, "> $outfile") || die "Can't create $outfile";
  print "creating $outfile\n";
  print F "$header\torientation\tpolyA\tpolyT\tpclipped\tnclipped\tprammate\tnrammate\teprammate\tenrammate\n";

  for my $s (keys %{$href}) {
    for my $chr (@{$href->{$s}{chr}}) {
      for my $site (@{$href->{$s}{$chr}{site}}) {
        my @pos = split(":", $site);
        my $hhref = $href->{$s}{$chr}{$site};
				&is_polyA2($hhref);

				my $orientation = "NA";
				if ($hhref->{polyA} eq "polyA" && $hhref->{polyT} eq "-") { $orientation = "+" } 
				if ($hhref->{polyT} eq "polyT" && $hhref->{polyA} eq "-") { $orientation = "-" } 
				
        #print F "$s\t$chr\t$pos[0]\t$pos[1]\t$hhref->{type}\t";
        print F "$hhref->{rec}\t$orientation\t";
        print F "$hhref->{polyA}\t$hhref->{polyT}\t";
        print F "$hhref->{pcontig}\t$hhref->{ncontig}\t";
		print F "$hhref->{prcontig}\t$hhref->{nrcontig}\t";
		print F "$hhref->{eprcontig}\t$hhref->{enrcontig}\n";
      }
    }
  }
}

sub write_contig
{
	my ($outprefix, $href, $header) = @_;
	my $outfile = $outprefix.".contig";
	
	open(F, "> $outfile") || die "Can't create $outfile";
	print "creating $outfile\n";
	#print F "sample\tchr\tram_start\tram_end\trep.repeat\tpos_clipped_contigs\tneg_clipped_contigs\tpos_clipped_singlets\tneg_clipped_singlets\tpos_ram_mate_contigs\tneg_ram_mate_contigs\tpos_ram_mate_singlets\tneg_ram_mate_singlets\n";
	print F "$header\tpos_polyA\tneg_polyA\tpos_spolyA\tneg_spolyA\tpos_clipped_contigs\tneg_clipped_contigs\tpos_clipped_singlets\tneg_clipped_singlets\tpos_ram_mate_contigs\tneg_ram_mate_contigs\tpos_ram_mate_singlets\tneg_ram_mate_singlets\n";

	for my $s (keys %{$href}) {
		for my $chr (keys %{$href->{$s}}) {
			for my $site (keys %{$href->{$s}->{$chr}}) {
				my @pos = split(":", $site);
				my $hhref = $href->{$s}->{$chr}->{$site};

				# put "-" if no sequence
				if (@{$hhref->{pcontig}} == 0) { push(@{$hhref->{pcontig}}, "-") }
				if (@{$hhref->{ncontig}} == 0) { push(@{$hhref->{ncontig}}, "-") }
				if (@{$hhref->{psinglet}} == 0) { push(@{$hhref->{psinglet}}, "-") }
				if (@{$hhref->{nsinglet}} == 0) { push(@{$hhref->{nsinglet}}, "-") }

				if (@{$hhref->{prcontig}} == 0) { push(@{$hhref->{prcontig}}, "-") }
				if (@{$hhref->{nrcontig}} == 0) { push(@{$hhref->{nrcontig}}, "-") }
				if (@{$hhref->{prsinglet}} == 0) { push(@{$hhref->{prsinglet}}, "-") }
				if (@{$hhref->{nrsinglet}} == 0) { push(@{$hhref->{nrsinglet}}, "-") }

				#print F "$s\t$chr\t$pos[0]\t$pos[1]\t$hhref->{type}\t";
				print F "$hhref->{rec}\t";
				print F "$hhref->{ppolyA}\t$hhref->{npolyA}\t$hhref->{pspolyA}\t$hhref->{nspolyA}\t";
				print F join(",", @{$hhref->{pcontig}})."\t".join(",", @{$hhref->{ncontig}})."\t";
				print F join(",", @{$hhref->{psinglet}})."\t".join(",", @{$hhref->{nsinglet}})."\t";
				print F join(",", @{$hhref->{prcontig}})."\t".join(",", @{$hhref->{nrcontig}})."\t";
				print F join(",", @{$hhref->{prsinglet}})."\t".join(",", @{$hhref->{nrsinglet}})."\n";
			}
		}
	}
}

sub make_contig2 
{
	my ($href, $odir, $tdir, $ronly, $conly, $eonly, $assembler, $verbose) = @_;

	my (@seq, @qual, @pcontig, @ncontig, @psinglet, @nsinglet, $cmd);
	for my $s (keys %{$href}) {
		for my $chr (@{$href->{$s}{chr}}) {
			for my $site (@{$href->{$s}{$chr}{site}}) {
				my $hhref = $href->{$s}->{$chr}->{$site};
				my $fprefix = "$tdir/$s.$chr.$site";
				$hhref->{pcontig} = "";
				$hhref->{ncontig} = "";
				$hhref->{prcontig} = "";
				$hhref->{nrcontig} = "";
				$hhref->{eprcontig} = "";
				$hhref->{enrcontig} = "";

				if (!$ronly && !$eonly) {
					# don't run cap3 for a single read
					if (@{$hhref->{pseq}} <= 1) { 
							$hhref->{pcontig} = pop(@{$hhref->{pseq}});
					} else {
						&write_fasta("$fprefix.pos.fa", \@{$hhref->{pseq}}, \@{$hhref->{pqual}}); 
						$cmd = $assembler." $fprefix.pos.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
						system($cmd);
		        		&parse_assembler_out2("$fprefix.pos.fa", $hhref);
					}

					if (@{$hhref->{nseq}} <= 1) { 
							$hhref->{ncontig} = pop(@{$hhref->{nseq}});
					} else {
						&write_fasta("$fprefix.neg.fa", \@{$hhref->{nseq}}, \@{$hhref->{nqual}}); 
						$cmd = $assembler." $fprefix.neg.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";  
						system($cmd);
   			     		&parse_assembler_out2("$fprefix.neg.fa", $hhref);
					}
				}

				if (!$conly && !$eonly) {
					if (@{$hhref->{prseq}} <= 1) { 
							$hhref->{prcontig} = pop(@{$hhref->{prseq}});
					} else {
		        		&write_fasta("$fprefix.rpos.fa", \@{$hhref->{prseq}}, \@{$hhref->{prqual}});
   			     		$cmd = $assembler." $fprefix.rpos.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
        				system($cmd);
		        		&parse_assembler_out2("$fprefix.rpos.fa", $hhref);
					}

					if (@{$hhref->{nrseq}} <= 1) {
						$hhref->{nrcontig} = pop(@{$hhref->{nrseq}});
					} else {
   			     		&write_fasta("$fprefix.rneg.fa", \@{$hhref->{nrseq}}, \@{$hhref->{nrqual}});
        				$cmd = $assembler." $fprefix.rneg.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
        				system($cmd);
        				&parse_assembler_out2("$fprefix.rneg.fa", $hhref);
					}
				}

			  if (!$conly && !$ronly) {
					print "eprseq: ".@{$hhref->{eprseq}}."\n";
					print "enrseq: ".@{$hhref->{enrseq}}."\n";
					if (@{$hhref->{eprseq}} > 0) {
					if (@{$hhref->{eprseq}} <= 1) {
                            $hhref->{eprcontig} = pop(@{$hhref->{eprseq}});
                    } else {
                        &write_fasta("$fprefix.erpos.fa", \@{$hhref->{eprseq}}, \@{$hhref->{eprqual}});
                        $cmd = $assembler." $fprefix.erpos.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
                        system($cmd);
                        &parse_assembler_out2("$fprefix.erpos.fa", $hhref);
                    }
					}

                    if (@{$hhref->{enrseq}} > 0) {
                    if (@{$hhref->{enrseq}} <= 1) {
                        $hhref->{enrcontig} = pop(@{$hhref->{enrseq}});
                    } else {
                        &write_fasta("$fprefix.erneg.fa", \@{$hhref->{enrseq}}, \@{$hhref->{enrqual}});
                        $cmd = $assembler." $fprefix.erneg.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
                        system($cmd);
                        &parse_assembler_out2("$fprefix.erneg.fa", $hhref);
                    }
					}
                }


        		$cmd = "rm $fprefix.*"; #print $cmd."\n";
		        if (!$verbose) { system($cmd); }

				if ($verbose) {
					if (!$ronly && !$eonly) {
					print "$s:$chr:$site:pseq:".@{$hhref->{pseq}}."\n";
					print "$s:$chr:$site:pqual:".@{$hhref->{pqual}}."\n";
					print "$s:$chr:$site:pcontig:".$hhref->{pcontig}."\n";

					print "$s:$chr:$site:nseq:".@{$hhref->{nseq}}."\n";
					print "$s:$chr:$site:nqual:".@{$hhref->{nqual}}."\n";
					print "$s:$chr:$site:ncontig:".$hhref->{ncontig}."\n";
					}
					if (!$conly && !$eonly) {
					print "$s:$chr:$site:prseq:".@{$hhref->{prseq}}."\n";
					print "$s:$chr:$site:prqual:".@{$hhref->{prqual}}."\n";
					print "$s:$chr:$site:prcontig:".$hhref->{prcontig}."\n";

					print "$s:$chr:$site:nrseq:".@{$hhref->{nrseq}}."\n";
					print "$s:$chr:$site:nrqual:".@{$hhref->{nrqual}}."\n";
					print "$s:$chr:$site:nrcontig:".$hhref->{nrcontig}."\n";
					}
					if (!$conly && !$ronly) {
					print "$s:$chr:$site:eprseq:".@{$hhref->{eprseq}}."\n";
					print "$s:$chr:$site:eprqual:".@{$hhref->{eprqual}}."\n";
					print "$s:$chr:$site:eprcontig:".$hhref->{eprcontig}."\n";

					print "$s:$chr:$site:enrseq:".@{$hhref->{enrseq}}."\n";
					print "$s:$chr:$site:enrqual:".@{$hhref->{enrqual}}."\n";
					print "$s:$chr:$site:enrcontig:".$hhref->{enrcontig}."\n";
					}
				}
			}
		}
	}
}


sub make_contig
{
	my ($href, $odir, $ronly, $conly, $assembler, $verbose) = @_;

	my (@seq, @qual, @pcontig, @ncontig, @psinglet, @nsinglet, $cmd);
	for my $s (keys %{$href}) {
		for my $chr (keys %{$href->{$s}}) {
			for my $site (keys %{$href->{$s}->{$chr}}) {
				my $hhref = $href->{$s}->{$chr}->{$site};
				my $fprefix = "$odir/$s.$chr.$site";
				@{$hhref->{pcontig}} = ();
				@{$hhref->{ncontig}} = ();
				@{$hhref->{psinglet}} = ();
				@{$hhref->{nsinglet}} = ();

				@{$hhref->{prcontig}} = ();
				@{$hhref->{nrcontig}} = ();
				@{$hhref->{prsinglet}} = ();
				@{$hhref->{nrsinglet}} = ();

				if (!$ronly) {
					# don't run cap3 for a single read
					if (@{$hhref->{pseq}} <= 1) { 
							@{$hhref->{pcontig}} = @{$hhref->{pseq}};
					} else {
						print "writing $fprefix.pos.fa\n";
						&write_fasta("$fprefix.pos.fa", \@{$hhref->{pseq}}, \@{$hhref->{pqual}}); 
						$cmd = $assembler." $fprefix.pos.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";  
						print $cmd."\n";
						system($cmd);
						print "parsing $fprefix.pos.fa \n";
        		&parse_assembler_out("$fprefix.pos.fa", $hhref);
						print "after parsing: pcontig ".join(",", @{$hhref->{pcontig}})."\n";
						print "after parsing: psinglet".join(",", @{$hhref->{psinglet}})."\n";
					}

					if (@{$hhref->{nseq}} <= 1) { 
							@{$hhref->{ncontig}} = @{$hhref->{nseq}};
					} else {
						print "writing $fprefix.neg.fa \n";
						&write_fasta("$fprefix.neg.fa", \@{$hhref->{nseq}}, \@{$hhref->{nqual}}); 
						$cmd = $assembler." $fprefix.neg.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";  
						print $cmd."\n";
						system($cmd);
        		&parse_assembler_out("$fprefix.neg.fa", $hhref);
						print "after parsing: ncontig ".join(",", @{$hhref->{ncontig}})."\n";
						print "after parsing: nsinglet ".join(",", @{$hhref->{nsinglet}})."\n";
					}
				}

				if (!$conly) {
					if (@{$hhref->{prseq}} <= 1) { 
							@{$hhref->{prcontig}} = @{$hhref->{prseq}};
					} else {
						print "writing $fprefix.rpos.fa \n";
        		&write_fasta("$fprefix.rpos.fa", \@{$hhref->{prseq}}, \@{$hhref->{prqual}});
        		$cmd = $assembler." $fprefix.rpos.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
						print $cmd."\n";
        		system($cmd);
        		&parse_assembler_out("$fprefix.rpos.fa", $hhref);
						print "after parsing prcontig".join(",", @{$hhref->{prcontig}})."\n";
						print "after parsing prsinglet".join(",", @{$hhref->{prsinglet}})."\n";
					}
					if (@{$hhref->{nrseq}} <= 1) { 
							@{$hhref->{nrcontig}} = @{$hhref->{nrseq}};
					} else {
						print "writing $fprefix.rneg.fa \n";
        		&write_fasta("$fprefix.rneg.fa", \@{$hhref->{nrseq}}, \@{$hhref->{nrqual}});
        		$cmd = $assembler." $fprefix.rneg.fa -i 21 -j 31 -o 16 -s 251 -p 70 > $fprefix.cap.log";
						print $cmd."\n";
        		system($cmd);
        		&parse_assembler_out("$fprefix.rneg.fa", $hhref);
						print "after parsing nrcontig".join(",", @{$hhref->{nrcontig}})."\n";
						print "after parsing nrsinglet".join(",", @{$hhref->{nrsinglet}})."\n";
					}
				}

				#if (!$ronly) 	{ &is_polyA($hhref) };

        $cmd = "rm $fprefix.*"; #print $cmd."\n";
        if (!$verbose) { system($cmd); }

				if ($verbose) {
					if (!$ronly) {
					print "$s:$chr:$site:pseq:".@{$hhref->{pseq}}."\n";
					print "$s:$chr:$site:pqual:".@{$hhref->{pqual}}."\n";
					print "$s:$chr:$site:pcontig:".@{$hhref->{pcontig}}."\n";
					print "$s:$chr:$site:psinglet:".@{$hhref->{psinglet}}."\n";

					print "$s:$chr:$site:nseq:".@{$hhref->{nseq}}."\n";
					print "$s:$chr:$site:nqual:".@{$hhref->{nqual}}."\n";
					print "$s:$chr:$site:ncontig:".@{$hhref->{ncontig}}."\n";
					print "$s:$chr:$site:nsinglet:".@{$hhref->{nsinglet}}."\n";
					}
if (!$conly) {
					print "$s:$chr:$site:prseq:".@{$hhref->{prseq}}."\n";
					print "$s:$chr:$site:prqual:".@{$hhref->{prqual}}."\n";
					print "$s:$chr:$site:prcontig:".@{$hhref->{prcontig}}."\n";
					print "$s:$chr:$site:prsinglet:".@{$hhref->{prsinglet}}."\n";

					print "$s:$chr:$site:nrseq:".@{$hhref->{nrseq}}."\n";
					print "$s:$chr:$site:nrqual:".@{$hhref->{nrqual}}."\n";
					print "$s:$chr:$site:nrcontig:".@{$hhref->{nrcontig}}."\n";
					print "$s:$chr:$site:nrsinglet:".@{$hhref->{nrsinglet}}."\n";
					}
				}
			}
		}
	}
}

sub longest
{
	my @a = @_;
	my $maxl = 0; my $maxs = "";
	foreach (@a) {
		if (length($_) > $maxl) { $maxl = length($_); $maxs = $_; }
	}
	return $maxs;
}

sub parse_assembler_out2
{
	my ($fa, $href) = @_;

	my (@contig, @singlet);
	&get_contig_seq("$fa.cap.contigs", \@contig);
	&get_contig_seq("$fa.cap.singlets", \@singlet);
	#print(join(", ", @contig));
	#print(join(", ", @singlet));

	if ($fa =~ m/\.pos.fa/) {
		$href->{pcontig} = &longest(@contig);
		if ($href->{pcontig} eq "") { $href->{pcontig} = &longest(@singlet); }
		if ($href->{pcontig} eq "") { $href->{pcontig} = "-"; }
	} elsif ($fa =~ m/\.neg.fa/) {
		$href->{ncontig} = &longest(@contig);
		if ($href->{ncontig} eq "") { $href->{ncontig} = &longest(@singlet); }
		if ($href->{ncontig} eq "") { $href->{ncontig} = "-"; }
	} elsif ($fa =~ m/\.rpos.fa/) {
		$href->{prcontig} = &longest(@contig);
		if ($href->{prcontig} eq "") { $href->{prcontig} = &longest(@singlet); }
		if ($href->{prcontig} eq "") { $href->{prcontig} = "-"; }
	} elsif ($fa =~ m/\.rneg.fa/) {
		$href->{nrcontig} = &longest(@contig);
		if ($href->{nrcontig} eq "") { $href->{nrcontig} = &longest(@singlet); }
		if ($href->{nrcontig} eq "") { $href->{nrcontig} = "-"; }
	} elsif ($fa =~ m/\.erpos.fa/) {
		$href->{eprcontig} = &longest(@contig);
		if ($href->{eprcontig} eq "") { $href->{eprcontig} = &longest(@singlet); }
		if ($href->{eprcontig} eq "") { $href->{eprcontig} = "-"; }
	} elsif ($fa =~ m/\.erneg.fa/) {
		$href->{enrcontig} = &longest(@contig);
		if ($href->{enrcontig} eq "") { $href->{enrcontig} = &longest(@singlet); }
		if ($href->{enrcontig} eq "") { $href->{enrcontig} = "-"; }
	}
}

sub parse_assembler_out
{
	my ($fa, $href) = @_;

	if ($fa =~ m/\.pos.fa/) {
		&get_contig_seq("$fa.cap.contigs", $href->{pcontig});
		&get_contig_seq("$fa.cap.singlets", $href->{psinglet});
		#print join(",", @{$href->{pcontig}})."\n";
		#print join(",", @{$href->{psinglet}})."\n";
	} elsif ($fa =~ m/\.neg.fa/) {
		&get_contig_seq("$fa.cap.contigs", $href->{ncontig});
		&get_contig_seq("$fa.cap.singlets", $href->{nsinglet});
		#print join(",", @{$href->{ncontig}})."\n";
		#print join(",", @{$href->{nsinglet}})."\n";
	} elsif ($fa =~ m/\.rpos.fa/) {
		&get_contig_seq("$fa.cap.contigs", $href->{prcontig});
		&get_contig_seq("$fa.cap.singlets", $href->{prsinglet});
		#print join(",", @{$href->{prcontig}})."\n";
		#print join(",", @{$href->{prsinglet}})."\n";
	} elsif ($fa =~ m/\.rneg.fa/) {
		&get_contig_seq("$fa.cap.contigs", $href->{nrcontig});
		&get_contig_seq("$fa.cap.singlets", $href->{nrsinglet});
		#print join(",", @{$href->{nrcontig}})."\n";
		#print join(",", @{$href->{nrsinglet}})."\n";
	}
}

sub get_contig_seq
{
	my ($infile, $ar) = @_;
	open(F, "< $infile") || die "Can't open $infile";
	my $l = "";
	while (<F>) { 
		chomp;
		if (m/^>/) {
			if ($l ne "") { push(@{$ar}, $l) }
			$l = "";
		} else {
			$l = $l.$_;
		}
	}
	if ($l ne "") { push(@{$ar}, $l) }
	close F;
}

sub is_polyA2 
{
  my ($href) = @_;

  my $seq = $href->{pcontig}; my ($acnt, $tcnt) = (0,0);
	if ($href->{pcontig} ne "-")  { $acnt =  &is_polyA_seq(substr($href->{pcontig}, -6, 6), "A") }
	if ($href->{ncontig} ne "-")  { $tcnt =  &is_polyA_seq(substr($href->{ncontig}, 0, 6), "T") }
	if ($acnt == 0) { $href->{polyA} = "-" } else { $href->{polyA} = "polyA" }
	if ($tcnt == 0) { $href->{polyT} = "-" } else { $href->{polyT} = "polyT" }
}

sub is_polyA
{
	my ($href) = @_;

	my $pcontig = ${$href->{pcontig}}[0];
	my $ncontig = ${$href->{ncontig}}[0];

	print "checking polyA for $pcontig : $ncontig\n";
	my ($pflag, $psflag, $nflag, $nsflag ) = ("-", "-", "-", "-"); # for +/-contigs and +/-singlets
	if (length($pcontig) >= 6) {
		
	if (&is_polyA_seq(substr($pcontig, 0, 6), "A") || &is_polyA_seq(substr($pcontig, -6, 6), "A")) {
		$pflag = "A";
	} elsif (&is_polyA_seq(substr($pcontig, 0, 6), "T") || &is_polyA_seq(substr($pcontig, -6, 6), "T") ) {
		$pflag = "T";
	}
	}

	if (length($ncontig) >=6) {
	if (&is_polyA_seq(substr($ncontig, 0, 6), "A") || &is_polyA_seq(substr($ncontig, -6, 6), "A")) {
		$nflag = "A";
	} elsif (&is_polyA_seq(substr($ncontig, 0, 6), "T") || &is_polyA_seq(substr($ncontig, -6, 6), "T") ) {
		$nflag = "T";
	}
	}

	foreach (@{$href->{psinglet}}) {
		my @a = split(/,/, $_);
		for (my $i=0; $i<@a; $i++) {
			if (length($a[$i]) >= 6) {
			 if (&is_polyA_seq(substr($a[$i], 0, 6), "A") || &is_polyA_seq(substr($a[$i], -6, 6), "A")) {
    			$psflag = "A";
  			} elsif (&is_polyA_seq(substr($pcontig, 0, 6), "T") || &is_polyA_seq(substr($pcontig, -6, 6), "T") ) {
    			$psflag = "T";
				}				
			}
  	}
	}

	foreach (@{$href->{nsinglet}}) {
		my @a = split(/,/, $_);
		for (my $i=0; $i<@a; $i++) {
			if (length($a[$i]) >= 6) {
			 if (&is_polyA_seq(substr($a[$i], 0, 6), "A") || &is_polyA_seq(substr($a[$i], -6, 6), "A")) {
    			$nsflag = "A";
  			} elsif (&is_polyA_seq(substr($pcontig, 0, 6), "T") || &is_polyA_seq(substr($pcontig, -6, 6), "T") ) {
    			$nsflag = "T";
				}				
			}
  	}
	}
	$href->{ppolyA} = $pflag;
	$href->{npolyA} = $nflag;
	$href->{pspolyA} = $psflag;
	$href->{nspolyA} = $nsflag;
}

sub is_polyA_seq
{
	my ($seq, $chr) = @_; 
	# A6 with at least five A
	my $cnt = () = ($seq =~ /$chr/g);
	
	if ($cnt >=5) { 
		return ($cnt) 
	} else {
		return(0)
	}
}

sub write_fasta
{
	my ($fprefix, $sref, $qref) = @_;  
	open(F, "> $fprefix") || die "Can't create $fprefix";
	for (my $i=0; $i< scalar @$sref; $i++) { print F ">cr$i\n$sref->[$i]\n"; }; 
	close F;

	#open(F, "> $fprefix.qual") || die "Can't create $fprefix.qual";
	#for (my $i=0; $i< scalar @$qref; $i++) { print F ">cr$i\n$qref->[$i]\n"; }; 
	#close F;
}

sub bam2fq
{
	my $usage = qq{
  Usage:   bam2fq [-p] [-v] <bam file> <out dir>

	};

	my %opts = ();
	getopts("hvp", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
	my $pair = 0; $pair = 1 if (defined($opts{p}));
	my $fname = shift(@ARGV);
	my $dir = shift(@ARGV);
	
	if (! -d $dir) { system("mkdir -p $dir") }
	
	my ( $outfn1, $outfn2, $outfn, $z1, $z2, $z, $out, %out_1, %out_2 );
	
	open(F, "$samtools view $fname |") or die "Can't open $fname: $!";
	my $prefix = basename($fname);


	my @str = split /\./, $prefix;


	$prefix = join(".", @str[0..($#str-1)]);
	
	#$prefix =~ s/.bam//;
	
	if ($verbose) { print "reading $fname ...\n";	}

	if ( $pair )
	{
	  $outfn1 = "$dir/$prefix\_1.fastq.gz";
	  $outfn2 = "$dir/$prefix\_2.fastq.gz";
	  $outfn = "$dir/$prefix.fastq.gz";
	
	  $z1 = new IO::Compress::Gzip $outfn1;
	  $z2 = new IO::Compress::Gzip $outfn2;
	  $z = new IO::Compress::Gzip $outfn;
	
		if ($verbose) { print "converting to \n $outfn1\n $outfn2\n $outfn ...\n";	}
	}
	else
	{
	  $outfn = "$dir/$prefix.fastq.gz";
	  $z = new IO::Compress::Gzip $outfn;

		if ($verbose) { print "converting to $outfn ...\n"; }	
	}
	
	while (<F>)
	{
	  chomp; my @F = split(/\t/);
	  if ($F[1] & 0x10)
	  {
	    ($F[9] = reverse $F[9]) =~ tr/gatcGATC/ctagCTAG/;
	    $F[10] = reverse $F[10];
	  }
	
	  $out = "\@$F[0]\n$F[9]\n+\n$F[10]\n";
	
	  if ( $pair )
	  {
        if ( ($F[1] & 1) && ($F[1] & 0x40) )
        { # first in pair
          if (exists $out_2{$F[0]})
          {    
            print $z1 $out;
            print $z2 $out_2{$F[0]};
            delete $out_2{$F[0]};
          }    
          else 
          {    
            $out_1{$F[0]} = $out;
          }    
        }    
        elsif ( ($F[1] & 1) && ($F[1] & 0x80) )
        {   # second in pair
          if (exists $out_1{$F[0]})
          {    
            print $z1 $out_1{$F[0]};
            print $z2 $out;
            delete $out_1{$F[0]};
          }    
          else 
          {    
            $out_2{$F[0]} = $out;
          }
        }    
        else 
        { # all the other reads
          print $z $out;
          if (exists $out_1{$F[0]})
          {    
            print $z $out_1{$F[0]};
            delete $out_1{$F[0]};
          }    
          elsif (exists $out_2{$F[0]})
          {    
            print $z $out_2{$F[0]};
            delete $out_2{$F[0]};
          }    
        }    
      }    
      else 
      {    
        print $z $out;
      }    
  }
 
  # print the orphans
  while ( my ($key, $val) = each %out_1 ) { print $z $val; }
  while ( my ($key, $val) = each %out_2 ) { print $z $val; }

  if ($verbose) { print "done bam2fq converting.\n"; }
}

sub get_chrl 
{
  my ($chrl_ref, $chrl_href, $organism) = @_;
  if ($organism eq "human" || $organism eq "hg18" || $organism eq "hg19") {
    #push(@$chrl_ref, "10"); $chrl_href->{10} = 1;
    for (my $i=1; $i<=22; $i++) { push(@$chrl_ref, "$i"); $chrl_href->{$i} = 1 };
    push(@$chrl_ref, "X"); $chrl_href->{"X"} = 1;
    push(@$chrl_ref, "Y"); $chrl_href->{"Y"} = 1;
  } elsif ($organism eq "orangutan" || $organism eq "ponAbe2") {
    for (my $i=1; $i<=22; $i++) {
      if ($i==2) { push(@$chrl_ref, "2a"); push(@$chrl_ref, "2b") }
      else { push(@$chrl_ref, "$i"); }
    }
    push(@$chrl_ref, "X");
    foreach my $c (@$chrl_ref) {
      $chrl_href->{$c} = 1
    }
  } elsif ($organism eq "chimpanzee" || $organism eq "panTro3") {
    for (my $i=1; $i<=22; $i++) {
      if ($i==2) { push(@$chrl_ref, "2A"); push(@$chrl_ref, "2B") }
      else { push(@$chrl_ref, "$i"); }
    }
    push(@$chrl_ref, "X");
    push(@$chrl_ref, "Y");
    foreach my $c (@$chrl_ref) {
      $chrl_href->{$c} = 1
    }
  } elsif ($organism eq "rhesus" || $organism eq "rheMac2") {
    for (my $i=1; $i<=20; $i++) { push(@$chrl_ref, "$i"); }
    push(@$chrl_ref, "X");
    foreach my $c (@$chrl_ref) {
      $chrl_href->{$c} = 1
    }
  }
}

sub cbam
{
	my $usage = qq{
  Usage:   cbam [options] <bam file>

  Options: -v     verboase 
           -c INT quality score cutoff (0: # from the bwa)
           -o STR outdir (.)
	};
	
	my %opts = ();
	getopts("hvc:o:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
	my $qcutoff = 0; $qcutoff = $opts{c} if (defined $opts{c});
	my $outdir = "."; $outdir = $opts{o} if (defined $opts{o});

	my $bamf = shift(@ARGV); 
	print "$bamf:\n";
	my $outf2 = basename($bamf);
	my $outf3 = basename($bamf);
	
	my @str = split /\./, $outf2;
	$outf2 = join(".", @str[0..($#str-1)]);
	$outf2 = $outf2 . '.softclips.consd.raw.bam';
	$outf3 = $outf2 . '.softclips.consd.cpos.bz2';
	
	#$outf2 =~ s/\.(sorted\.|)bam/.softclips.consd.raw.bam/;
	$outf2 = "$outdir/$outf2";  
	#$outf3 =~ s/\.(sorted\.|)bam/.softclips.consd.cpos.bz2/;
	$outf3 = "$outdir/$outf3";  

	open(O2, "| $samtools view -bS - > $outf2") || die "Can't create $outf2";
	open(O3, "| bzip2 - > $outf3") || die "Can't create $outf3";

	# print the bam header
  open(F, "$samtools view -H $bamf |");
	my @header = <F>;
	print O2 @header;
  close F;

	open(F, "$samtools view $bamf|") || die "Can't open $bamf";
	my ($cigar, $qual, $qual1, $qual2, $s1, $s2, $selected, $selected2, $cpos);
	while(<F>) {
		chomp; my @a=split(/\t/);
    if ($a[1] & 0x0004) { next } # unmapped ; there must not be this case for the clipped reads, but just in case
		$cigar = $a[5]; $qual = $a[10];
		$s1 = $cigar;
		$s2 = $cigar;
		$selected = $selected2= 0;	
    if ($s1 =~ /^(\d+)S\S+/) { # clippeding in the beginning 
        if ($1 >= 5) {
        	$qual1 = substr($qual, 0, $1);
        	if ($qual1 =~ m/([^#]+)/) {
						 if (length($1) >= 5) {  # more than 5 good quality bases 
							$selected = 1;	
							if (! ($a[1] & 0x0010)) { # positive strand mapping 
								$selected2 = 1; 
								$cpos = $a[3];
							}
						}
					}
				}
    }
    if ($s2=~ /\D(\d+)S$/) { # clipping in the end
      if ($1 >= 5) {
 	    	$qual2 = substr($qual, -$1);
 	     	if ($qual2 =~ m/([^#]+)/) { 
					if (length($1) >= 5 ) { 
						$selected = 1;	
						if ($a[1] & 0x0010) { # negative strand mapping
							$selected2 = 1; 
						  $cpos = &get_cpos($a[3], $cigar, $qual2, -1);
						} 
					}
				} 
			}
    }
		if ($selected2) {
				print O2 "$_\n"; 
				print O3 "$a[2]\t$cpos\n";
		}
	}
	close F;
	close O2;
	close O3;

	print "done generating cbam: $outf2\n";

	print "start sorting and generating the index for $outf2 ...\n";
	my $prefix = $outf2;
	$prefix =~ s/\.raw\.bam$//;
	system("$samtools sort $outf2 $prefix; $samtools index $prefix.bam; rm $outf2");
	print "done generating cbam and its index.\n";
	
}

sub get_cpos
{
	my ($pos, $cigar, $qual, $strand) = @_;
	my $cpos = $pos;
	#print "pos: $pos, cigar: $cigar, qual: $qual, strand: $strand\n";
		
	if ($strand == 1) { return($cpos) }

	# for the negative strand clipped reads, count the number of base pairs (M & I not D)
	my $gap = 0;
	$cigar =~ s/^\d+S(.*)/$1/;
	my @cnt = split(/[MINHPSD]/, $cigar);
	my @chr = split(/\d+/, $cigar);
	shift @chr;
	for (my $i=0; $i<@chr; $i++) {
		if ($chr[$i] ne "I" && $chr[$i] ne "S") { $gap += $cnt[$i] }
	}
	$cpos = $cpos + $gap;
	#print "gap:$gap\n";
	#print "cpos: $cpos\n";	
	return(-$cpos);
}

sub cisize
{
	my $usage = qq{
  Usage:   cisize <input isize file>
	
	};

	if (@ARGV < 1) { die $usage }
	my $file = shift(@ARGV);
	open(F, "< $file") || die "can't open $file";
	
	my $outf1 = $file; # eg. gbm0145_cancer.isinfo
	$outf1 =~ s/\.isinfo/\.isize/;
	my $outf2 = $file; 
	$outf2 =~ s/\.isinfo/\.rl/; 

	open(O1, "> $outf1") || die "can't create file $outf1";
	open(O2, "> $outf2") || die "can't create file $outf2";

	my (%is_mean, %is_median, %is_sd, %rl);
	my ($rg, $val);
	while (<F>) {
		chomp;
		if (m/Read length/) {	
			$rg =$_;
			$rg =~ s/Read length:\t(\S+)$/$1/; 
			$val = <F>; chomp($val); 
			$rl{$rg} = $val;
			#print "wrote $val into rl hash with key $rg\n";
		} elsif (m/Mean insert size/) {
			$rg = $_;
			$rg =~ s/Mean insert size:\t(\S+)$/$1/; 
			$val = <F>; chomp($val);
			$is_mean{$rg} = $val;
			#print "wrote $val into is_mean hash with key $rg\n";
		} elsif (m/Median insert size/) {
			$rg = $_;
			$rg =~ s/Median insert size:\t(\S+)$/$1/; 
			$val = <F>; chomp($val);
			$is_median{$rg} = $val;
			#print "wrote $val into is_median hash with key $rg\n";
		} elsif (m/Standard deviation of insert size/) {
			$rg = $_;
			$rg =~ s/Standard deviation of insert size:\t(\S+)$/$1/; 
			$val = <F>; chomp($val); 
			$is_sd{$rg} = $val;
			#print "wrote $val into is_sd hash with key $rg\n";
		}
	}

	# representative values (mean) for all read groups
	my $all_is_mean = int(&mean(&get_hash_values(\%is_mean)));
	my $all_is_median = int(&mean(&get_hash_values(\%is_median)));
	my $all_is_sd = int(&mean(&get_hash_values(\%is_sd)));
	my $all_rl = int(&mean(&get_hash_values(\%rl)));

	for my $rg (keys %is_mean) { print O1 "rg$rg\t".int($is_mean{$rg})."\t".int($is_sd{$rg})."\n"; } 
	print O1 "all\t$all_is_mean\t$all_is_sd\n";

	while (my ($key, $val) = each %rl) { print O2 "rg$key\t$val\n"; }
	print O2 "all\t$all_rl\n";
}

sub get_hash_values
{
	my $href = @_[0];
	my @val = ();
	while (my ($key, $val) = each (%$href)) {
		push(@val, $val);
	}
	return (\@val);
}

sub mean
{
	@_ == 1 or die ('Sub usage: $average = average(\@array);'); 
	my ($array_ref) = @_; 
	my $sum; 
	my $count = scalar @$array_ref; 
	foreach (@$array_ref) { $sum += $_; } 
	return $sum / $count; 
}

sub median 
{
	@_ == 1 or die ('Sub usage: $median = median(\@array);');
	my ($array_ref) = @_;
	my $count = scalar @$array_ref;
	# Sort a COPY of the array, leaving the original untouched
	my @array = sort { $a <=> $b } @$array_ref;
	if ($count % 2) {
	return $array[int($count/2)];
	} else {
		return ($array[$count/2] + $array[$count/2 - 1]) / 2;
	}
}

sub bamram
{
  my $usage = qq{
  Usage:   bamram [options] <ref bam> <rambam> <ram bz2> <rbam1> <rbam2>

  Options: -v         verbose
           -d STR     out dir

  };

	my %opts = ();
	getopts("hvd:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
	my $dir = "."; $dir = $opts{d} if (defined $opts{d});

	my $refbam = "$dir/".shift(@ARGV);
	my $rbamf = "$dir/".shift(@ARGV);
	my $ramf = "$dir/".shift(@ARGV);

	open(REF, "$samtools view $refbam | ") or die "no reference bam: $refbam: $!";

	# load repeat mappings into a hash
	my %h = (); # record end(1/2):READNAME\tREPAT_NAME1:fr1,REPEAT_NAME2:fr2, ...
	my ($base, $rasym); 
	my $exo=0; 	

	while (@ARGV) {
		my $rabam = "$dir/".shift(@ARGV);
		if ($rabam =~ m/\.um/) { 
			$exo = 1; # exogeneous
			print "exogenous ..\n";
		}

		my $end = $rabam;
		$end =~ s/.*_(1|2).*\.bam/$1/;
		print "reading bamfile: $rabam, end: $end\n";
		my ($mapcnt, $dupmapcnt) = (0, 0);
		open(RA, "$samtools view $rabam | ") or die "no reapeat assembly  bam: $rabam: $!";
		while (<RA>) {
				chomp; my @a=split(/\t/);
				if ($a[1] & 0x4) { next; }
				my %c = ();
				my ($key, $value);
			 	if (m/XA:Z:(\S+)/) {
					my @b = split(/;/, $1); 
					for (my $i=0; $i<@b; $i++) {
						my @rnames = split(/,/, $b[$i]);
						if (exists $c{$rnames[0]}) { $c{$rnames[0]}++; }
						else { $c{$rnames[0]} = 1 }
					}
				}
				$key = "$end:$a[0]";	
				$value = $a[2];	
				my @rmap = ();
				while (my ($k, $v) = each %c) {
					push(@rmap,"$k:$v"); 
				}
				if (@rmap >= 2) {
					$dupmapcnt++;
					$value = $value.",".join(",", @rmap);
				}
				
				$h{$key} = $value;  
				#print "$key\t$value\n";
				$mapcnt++;
		}
		print "done processing $rabam: mapcnt: $mapcnt: dupmapcnt: $dupmapcnt\n";
	}

	open(O, "| bzip2 - > $ramf") or die "can't create $ramf";
	open(O2, "| $samtools view -bS - > $rbamf") or die "can't create $rbamf";

  # print the bam header
	print "starting reading $refbam\n";
  open(F, "$samtools view -H $refbam |");
  my @header = <F>;
	print O2 @header;
  close F;

	# matching read ids from the referent bam and generate a ram file and a bam only with rams
	open(F, "$samtools view $refbam |" )or die "can't open $refbam";
	my ($cnt, $str, $map, $prname) = (0, "", "", "");

	my %ram = (); # in the case of cl=1, generate unique RAM map positions out of max two read pairs
	              # record readname => chr:pos:rname to generate unique pos 
                # for the two pairs of reads originated from one pair
                # e.g. HWI-ST1001:7:2210:14748:43848mu1 and HWI-ST1001:7:2210:14748:43848mu2
                # with the exactly same mpos to the same repeat type

	print "starting the id matching to generate a bam from a refbam ...\n";

	while (<F>) {
		chomp; 
		my @a = split(/\t/);
		$a[2] =~ s/chr//; # drop chr
    	#if (!exists($chrl{$a[2]}) || $a[1] & 0x4 || (m/XT:A:/ && !m/XT:A:U/) || $a[4] == 0 ) { next; }
    	if ($a[1] & 0x4 || (m/XT:A:/ && !m/XT:A:U/) || $a[4] == 0 ) { next; }
    	if ($a[1] & 0x40) {# end:1 
			if (exists $h{"2:$a[0]"}) {
  	   			 if ($a[1] & 0x10) { 
					$map = "$a[2]\t-$a[3]\t".$h{"2:$a[0]"}
				} else {
					$map = "$a[2]\t$a[3]\t".$h{"2:$a[0]"}
				}
				# check whether the record with the same readname and ram pos exists 
				# (only differ in the last suffix 1 or 2 for mu & sc)same readname exists
				$prname = $a[0];
				if ($prname =~ /mu1$/) { $prname =~ s/mu1/mu2/ }
				elsif ($prname =~ /mu2$/) { $prname =~ s/mu2/mu1/ }

				if ($exo == 0 || ($exo == 1 && ! exists $ram{$prname}) || ($exo==1 && $ram{$prname} ne $map)) {
					#debug begin
		         	#print "$a[0]\t$map\n";
					#debug end
		         	print O "$a[0]\t$map\n";
					if ($exo==1) { $ram{$a[0]} = $map }
					$cnt++;
					#print "end1:identified mate map:$_\n";
					$str = "$a[0]:\"".$h{"2:$a[0]"}."\"";
					$str = substr($str, 0, min(length($str), 250)); # trim rname when the te name is too long
					$_ =~ s/$a[0]/$str/;
					#debug begin
		         	#print "$_\n";
					#debug end
					print O2 "$_\n";
					#print "$cnt: $map to $ramf\n";
				}
			}
		} elsif ($a[1] & 0x80) { # end:2
			if (exists $h{"1:$a[0]"}) {
      			if ($a[1] & 0x10) { 
					$map = "$a[2]\t-$a[3]\t".$h{"1:$a[0]"}		
				} else {
					$map = "$a[2]\t$a[3]\t".$h{"1:$a[0]"}		
				}
				$prname = $a[0];
				if ($prname =~ /mu1$/) { $prname =~ s/mu1/mu2/ }
				elsif ($prname =~ /mu2$/) { $prname =~ s/mu2/mu1/ }
				if ($exo == 0 || ($exo == 1 && ! exists $ram{$prname}) || ($exo==1 && $ram{$prname} ne $map)) {
					#debug begin
		         	#print "$a[0]\t$map\n";
					#debug end
	   		      	print O "$a[0]\t$map\n";
					if ($exo==1) { $ram{$a[0]} = $map }
					$cnt++;
					#print "end2:identified mate map:$_\n";
					$str = "$a[0]:\"".$h{"1:$a[0]"}."\"";
					$str = substr($str, 0, min(length($str), 250)); # trim rname when the te name is too long
					$_ =~ s/$a[0]/$str/;
					#debug begin
		         	#print "$_\n";
					#debug end
					print O2 "$_\n";
					#print "$cnt: $map to $ramf\n";
				}
      		}
    	}
  	}
	close F;
	close O; 
	close O2;

	print "done generating ram and ram.bam with ".$cnt. " rams.\n";
}

sub ram 
{
	my $usage = qq{
  Usage:   ram [options] <ref bam file/dir> <repeat read names dir> <out dir>

  Options: -v     verboase 
           -d     input is folder 
           -s     start sram without running ram
           -q STR queue name (all_2h)

  Example: ram -d /files/CBMI/parklab/alee/ra/gbm/bam/06-0185-01A /files/CBMI/parklab/alee/ra/gbm/rname/06-0185-01A /files/CBMI/parklab/alee/ra/gbm/rami/06-0185-01A
	};

	my %opts = ();
	getopts("hvdsq:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $queue = "all_2h"; $queue = $opts{q} if (defined $opts{q});
  my $sram = 0; $sram = $opts{s} if (defined $opts{s});

	my $in_ref = shift(@ARGV);
	my $in_repeat = shift(@ARGV);
	my $out_dir = shift(@ARGV);

	my $rfile1 = "/groups/park/alee/ra/data/rmasker/hg18.repeats.txt.00"; # repeat names1
	my $rfile2 = "/groups/park/alee/ra/data/rmasker/hg18.repeats.txt.01"; # repeat names2
	 
  unless (-e $out_dir) { system("mkdir -p $out_dir") };

	if (defined $opts{d}) {
	unless ($sram) {
		my @files = <$in_ref/*.rg*.bam>;
		my %jobs = (); 
		for my $f (@files) {
      my $base = basename($f);
			$base =~ s/.bam//;
			my $jid = `bjobs -w | awk '/ram.$base/ { print \$1 }'`;
      if (-e "$out_dir/$base.ram.ok" || $jid ne "" ) {  next; }

     	my $cmd = "bsub -R \"rusage[mem=8000]\" -q $queue -J ram.$base -o $log_dir/ram.$base.out \"ra.pl ram -v $f $in_repeat $out_dir \" | perl -ne 'm/<(\\d+)>/; print \$1'";
       print "$cmd \n" if ($verbose);
       $jid = `$cmd`;
			 $jobs{$base} = $jid;	
		}
		my @jobs = keys(%jobs);
		print "\n\n".scalar(@jobs)." jobs submitted\n";

		my @completed = ();
		print "checking job completion\n" if ($verbose);
		while (1) {
      @completed = ();
      for my $f (@files) {
        my $base = basename($f); $base =~ s/.bam//;
        if (-e "$out_dir/$base.ram.ok") { push(@completed, $base); }
      }
      print scalar(@completed)." jobs completed\n" if ($verbose);
      if (scalar(@completed) < scalar(@files)) {
        sleep(60);
      } else {
        last;
      };
    }
		}

   print "checking sram(ram file per repeat)\n" if ($verbose);
		unless (-e "$out_dir/sram.".basename($rfile1).".ok") {
			my $jid1 = &bsub_sram($out_dir, $rfile1); # submit the job and write sram.ok when it's done
		}
		unless (-e "$out_dir/sram.".basename($rfile2).".ok") { 
			my $jid2 = &bsub_sram($out_dir, $rfile2);
		}

		my $completed=0;
   	while (1) {
    	print "checking sram job completion\n" if ($verbose);
			if ($completed >= 2) { 
				print "done generating sram files\n";
				last; 
			}
			if (-e "$out_dir/sram.".basename($rfile1).".ok") {
				$completed++;
				print "sram for rfile1  completed: $completed\n" if ($verbose);
			}
			if (-e "$out_dir/sram.".basename($rfile2).".ok") {
				$completed++;
				print "sram for rfile2  completed: $completed\n" if ($verbose);
			}
			sleep(60);
		}
  } else {
    &ram_bam($in_ref, $in_repeat, $out_dir, $verbose);
  }
}

sub bsub_sram 
{
	my ($out_dir, $rfile) = @_;
	my $queue = "park_12h";
	my $job = "sram.".basename($out_dir).".".basename($rfile);
	my $cmd = "bsub -q $queue -J $job -o $log_dir/$job.out \"ra.pl sram $out_dir $rfile\" | perl -ne 'm/<(\\d+)>/; print \$1'";
	my $jid = `$cmd`;
	print "$jid: $cmd\n";

	return($jid);
}

sub sram
{
	my $usage = qq{
  Usage:   sram [options] <source ram dir> <repeat list file> 

  Options: -v     verboase 
           -q STR queue name (park_12h)

	};

	my %opts = ();
	getopts("hvq:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $queue = "park_12h"; $queue = $opts{q} if (defined $opts{q});

	# out dir is the same as source ram dir
	my $out_dir = shift(@ARGV);
	my $rfile = shift(@ARGV);

	my %rpr = (); # ram file handles 

	#create file handle for each repeat
	open(R, "< $rfile") or die "can't open $rfile";
	my @r = <R>;
	for my $r (@r) {
		my @a = split(/\t/, $r);
    my $filename = "$out_dir/$a[0].sram.bz2";
		#print "creating $filename\n" if ($verbose);
    my $fh = new IO::File "| bzip2 - > \"$filename\"";
    $rpr{$a[0]} = $fh; 
  }
	print "done creating sram files.. \n" if ($verbose);

	#separte ram positions per repeat
	my @files = <$out_dir/TCGA*.ram.bz2>;
	my $fh; 
	print "reading ram files ..\n" if ($verbose);
	#for (my $i=0; $i<1; $i++) {
	for my $f (@files) {
		#my $f = $files[$i];
		open(F, "bzcat $f |") or die "can't open $f";
		print "reading ".basename($f)."\n" if ($verbose);
		while (<F>) {
			chomp;
			my @a = split(/\t/);
			if (exists $rpr{$a[0]}) {
				$fh = $rpr{$a[0]}; 
				if ($a[2] eq "+") { print $fh "$a[1]\t$a[3]\t$a[4]\n"; }
				else { print $fh "$a[1]\t$a[3]\t-$a[4]\n"; }
				#print $fh join("\t", @a));
			}
		}
		close F;
	}
	&close_fh(\%rpr);
	&rm_empty_sram(\@r, $out_dir); # remove empty sram files
	system("echo @files > $out_dir/sram.".basename($rfile).".ok");
}

sub close_fh
{
	my $fh_href = @_[0];
	for my $key (keys %$fh_href) {
		close $fh_href->{$key};
	}
}

sub rm_empty_sram
{
	my ($rref, $out_dir) = @_;
	for my $r (@$rref) {
		my @a = split(/\t/, $r);
		my $size = -s "$out_dir/$a[0].sram.bz2";
		if ($size == 0 || $size == 14) {
			system("rm \"$out_dir/$a[0].sram.bz2\"");
			print "done removing $out_dir/$a[0].sram.bz2\n";
		}
	}
}

sub fill_hash
{
	my ($f, $href, $end, $verbose) = @_;
	open(F, "bzcat $f |") or die "Can't open $f";
	print "reading ".basename($f)."\n" if ($verbose);
	while (<F>){
		chomp; my @a = split(/\t/);
		if (exists $href->{$a[0]}) {
			print "$a[0]:$a[1] already in the hash!: this shouldn't happen\n"
		}
		$href->{"$end:$a[0]"} = $a[1];	
	}
}


# reading through a bam and write the map pos for the ra reads
# print records to $out_dir/$id
sub match_id 
{
	# extract rg id and append it as a prefix for the map files ID#$a[0] 
	my ($bamf, $rg, $href, $out_dir, $verbose) = @_;
	print "matching ids for $bamf\n" if ($verbose);
	my $base = basename($bamf); $base =~ s/.bam//;
#	my $outf = "$out_dir/$base.ram.bz2";
	my $rbam = "$out_dir/$base.ram.bam";
#	open(O, "|bzip2 - > $outf") || die "Can't create $outf";

	open(O2, "| $samtools view -bS - > $rbam") || die "Can't create $rbam";

	# print the bam header
  open(F, "$samtools view -H $bamf |");
	my @header = <F>;
	print O2 @header;
  close F;

	open(F, "$samtools view $bamf|") || die "Can't open $bamf";
	my $cnt=0;
	while (<F>) {
		chomp;
		my @a = split(/\t/);
		if ($a[1] & 0x4 || !(m/XT:A:U/) || $a[2] eq "M" || $a[2] eq "chrM") { next; }
		if ($a[1] & 0x40 ) {# end:1	
			if (exists $href->{"2:$a[0]"}) {
				print O2 "$_\n";
#				if ($a[1] & 0x10) { 
#					print O $href->{"2:$a[0]"}."\t$rg:$a[0]\t-\t$a[2]\t$a[3]\n";
#				} else {
#					print O $href->{"2:$a[0]"}."\t$rg:$a[0]\t+\t$a[2]\t$a[3]\n";
#				}
				$cnt++;
			}			
		} elsif ($a[1] & 0x80) { # end:2
			if (exists $href->{"1:$a[0]"}) {
				print O2 "$_\n";
#				if ($a[1] & 0x10) { 
#					print O $href->{"1:$a[0]"}."\t$rg:$a[0]\t-\t$a[2]\t$a[3]\n";
#				} else {
#					print O $href->{"1:$a[0]"}."\t$rg:$a[0]\t+\t$a[2]\t$a[3]\n";
#				}
				$cnt++;
			}			
		}
	}
	close F;
#	close O;
	system("echo $cnt > $out_dir/$base.ram.ok"); # will be map cnts 
}

sub rbam_bam
{
	my ($bamf, $in_repeat, $out_dir, $verbose) = @_;
	my %h = ();

	my $base = basename($bamf);
	$base =~ s/.bam//;
	my $id = $base;
	$id =~ s/\.(rg.+)//;
	my $rg = $1;
	print "$base: $id: $rg\n" if ($verbose);

	my @rfiles = <$in_repeat/${id}_?.$rg.fastq*.names.bz2>; 
	#my @rfiles = <$in_repeat/${id}.$rg.*.fastq*.names.bz2>; 
	print "in_repeat: $in_repeat, id: $id, rg: $rg, rfiles: @rfiles\n";
	for my $f (@rfiles) {
		$f =~ /${id}_(1|2).*/;
		my $end = $1; 
		&fill_hash($f, \%h, $end, $verbose); 
		#print keys(%h);
	} 
	# read bam and print repeat rname +/0 pos outfilename $base.ram 
	&make_rbam($bamf, $rg, \%h, $out_dir, $verbose);
}

sub make_rbam
{
	my ($bamf, $rg, $href, $out_dir, $verbose) = @_;
	print "matching ids for $bamf\n" if ($verbose);
	my $base = basename($bamf); $base =~ s/.bam//;
	my $rbam = "$out_dir/$base.ram.bam";
	open(O2, "| $samtools view -bS - > $rbam") || die "Can't create $rbam";

	# print the bam header
  open(F, "$samtools view -H $bamf |");
	my @header = <F>;
	print O2 @header;
  close F;

	open(F, "$samtools view $bamf|") || die "Can't open $bamf";
	my $cnt=0;
	while (<F>) {
		chomp;
		my @a = split(/\t/);
		if ($a[1] & 0x4 || !(m/XT:A:U/) || $a[2] eq "M" || $a[2] eq "chrM") { next; }
		if ($a[1] & 0x40 ) {# end:1	
			if (exists $href->{"2:$a[0]"}) {
				print O2 "$_\n";
				$cnt++;
			}			
		} elsif ($a[1] & 0x80) { # end:2
			if (exists $href->{"1:$a[0]"}) {
				print O2 "$_\n";
				$cnt++;
			}			
		}
	}
	close F;
	close O;
	system("echo $cnt > $out_dir/$base.ram.ok"); # will be map cnts 
}

sub ram_bam
{
	my ($bamf, $in_repeat, $out_dir, $verbose) = @_;
	my %h = ();

	my $base = basename($bamf);
	$base =~ s/.bam//;
	my $id = $base;
	$id =~ s/\.(rg.+)//;
	my $rg = $1;
	print "$base: $id: $rg\n" if ($verbose);

	#my @rfiles = <$in_repeat/${id}_?.$rg.fastq*.names.bz2>; 
	my @rfiles = <$in_repeat/${id}_?.$rg.*.fastq*.names.bz2>; 
	print "in_repeat: $in_repeat, id: $id, rg: $rg, rfiles: @rfiles\n";
	for my $f (@rfiles) {
		$f =~ /${id}_(1|2).*/;
		my $end = $1; 
		&fill_hash($f, \%h, $end, $verbose); 
		#print keys(%h);
	} 
	# read bam and print repeat rname +/0 pos outfilename $base.ram 
	&match_id($bamf, $rg, \%h, $out_dir, $verbose);
}

sub mapcnt
{
 my $usage = qq{
  Usage:    mapcnt [options] <sample id>
      
  Options: -d STR  base dir (.)
           -c      map soft-clipped reads (yes) 
           -p INT  number of threads for bwa alignment          
           -l INT  bwa aln seed lenegh (40)
           -k INT  bwa aln maximum differences in the seed (2)
           -n INT  bwa aln max #diff (int) or missing prob under 0.02 err rate (floot) (3)
           -r STR  sequence assembly (.fasta)
           -R STR  assembly symbol (ra and va for endogenoeus and exogenous)
           -s STR  scratch or temp space to save the intermediate mapping files (/scratch/el114)

    };

  my %opts = ();
    getopts("hd:cp:l:k:n:r:R:s:", \%opts);
  if (@ARGV < 1 || exists($opts{h})) { die $usage };

    my $sample = shift(@ARGV);

    # set defaults
	my $dir= ".";
    my $l = 40;
    my $k = 2;
    my $n = 3;

	my $fq = "$dir/$sample/bam/$sample.unmapped.fq.gz";
	my $scdir= "$dir/$sample/bam/$sample";
	my $teapl = "~/repeat_anlaysis/code/Tea/scripts/tea.pl";
	
	my ($queue, $job, $job_log, $assembly, $sai, $errfile, $umf, $outprefix);

	`bsub -q $queue -J $job -o $job_log "$bwa aln $assembly $fq > $sai 2>>$errfile;
     $bwa samse -n 100 $assembly $sai $fq 2>>$errfile | $teapl vrs  -o $outprefix -f $umf;
	rm $sai"
	`;
}

sub vrs
{
    my %opts = ();
    getopts("o:f:", \%opts);
    die("Usage: vrs.pl [-o outprefix -f fastq] \n") unless (defined($opts{o})  && defined($opts{f}));

    my $outprefix = $opts{o};
    my $umf = $opts{f};

#   if (! -d "$outprefix/bam") {    system("mkdir -p $outprefix/bam") }; 
#   if (! -d "$outprefix/rnames") { system("mkdir -p $outprefix/rnames") }; 
    if (! -d "$outprefix/rcnt") {   system("mkdir -p $outprefix/rcnt") };

    my $rnamef = $umf;
    my $rcntf = $umf;
    my $bamf = $umf;
    $rnamef =~ s/fq.gz/rnames.gz/;
    $rcntf =~ s/fq.gz/rcnt/;
    $bamf =~ s/fq.gz/bam/;
    $rnamef = "$outprefix/rnames/$rnamef";
    $rcntf = "$outprefix/rcnt/$rcntf";
    $bamf = "$outprefix/bam/$bamf";

#   open(RN, "| gzip -c > $rnamef") || die "can't create $rnamef";
    open(RC, "> $rcntf") || die "can't create $rcntf";
#   open(B, "| samtools view -bS - > $bamf") || die "can't create $bamf";

    my @seqs = (); # the list of seq libraries 
    my %rnames = (); # virus name -> an array to read names 
    while (<>) {
        if (m/^\@SQ/) {
            m/SN:(\S+)/;
            my $sn = $1;
            push(@seqs, $sn);
        }
#       print B $_;
        chomp;
        my @a = split(/\t/);
        if ($a[2] ne "*") { # mapped to the virus assembly
            my ($rname, $vname) = ($a[0], $a[2]);
            $rname =~ s/_(1|2)$//; # trim the suffix

            if (m/XT:A:U/) {
                #print "added: $rname to $vname\n";
                push(@{$rnames{$vname}}, $rname);
            } elsif (m/XA:Z:(\S+)/) {
                my (@tokens, $altvname); # if a read is mapped to only one virus 
                my %alt = ();     # virus -> mapping count for this read 
                $alt{$vname} = 1;
                my @b = split(/;/, $1);
                for (my $i=0; $i<@b; $i++) {
                    my @tokens = split(/,/, $b[$i]);
                    my $altvname = $tokens[0];
                    if (exists $alt{$altvname}) { $alt{$altvname}++; }
                    else { $alt{$altvname} = 1 }
                }
                if ((my @keys = keys %alt) == 1) {
                    #print "added: $rname to $vname\n";
                    push(@{$rnames{$vname}}, $rname)
                }
            }
        }
    }

    # print the read names and counts
    for my $v (@seqs) {
        if (exists $rnames{$v}) {
#           print RN "$v\t@{$rnames{$v}}";
            print RC "$v\t".scalar @{$rnames{$v}}."\n";

        } else {
            print RC "$v\t0\n";
        }
    }
}

sub rmdup
{
    my ($sbamf, $obamf, $verbose) = @ARGV;
	$verbose = 0 if ($verbose eq "" );

    open(F, "$samtools view -hX $sbamf |") or die "can't open $sbamf";
    open(O, "| $samtools view -bS - > $obamf |") or die "can't create $obamf";
    my %h = (); # location -> uq and um read id and seq

    my (@a, $pos, $oldpos) = ((), "", "");
    while (<F>) {
        if (m/^@/) { print O; next }; # print header
        chomp; @a=split(/\t/);
        if ($a[1] =~ m/u/ && $a[1] =~ m/U/) { next }
        $pos = "$a[2]:$a[3]";
        my @x = keys %h;
        #print "$a[0]:$pos:hash: ".scalar @x."\n";
        if (exists $h{$pos}) { # the same locus appeared before 
            print "happend before\n" if $verbose;
            if ($a[1] =~ m/u/) { # unmapped read
                # check for the identical unmapped sequences 
                print "unmapped.. " if $verbose;
                my $dup = 0;
                for my $r (keys %{$h{$pos}}) {
                    if (exists $h{$pos}{$r}{urec}) {
                        my @b = split(/\t/, $h{$pos}{$r}{urec});
                        if ($a[9] eq $b[9]) {
                            $dup = 1;
                        }
                    }
                }
                if ($dup == 0) {
                    $h{$pos}{$a[0]}{urec} = $_;
                    print "no duplicate.. inseted into hash for locus $pos:$a[0] urec..\n" if $verbose;
                }
            } elsif ($a[1] =~ m/U/ && $a[4]>0 && (m/XT:A:U/ || !m/XT:A/)) {
                $h{$pos}{$a[0]}{mrec} = $_;
                print "mapped.. inserted into hash for $pos:$a[0] mrec..\n" if $verbose;
            }

        } else { # the new locus
            print "new locus.." if $verbose;
            print "oldpos: $oldpos\n" if $verbose;

            if (exists $h{$oldpos}) {
                print "hash exists\n" if $verbose;
                # print each pair if both mapped and unmapped read exist in the hash    
                for my $r (keys %{$h{$oldpos}}) {
                    if (exists $h{$oldpos}{$r}{mrec} && exists $h{$oldpos}{$r}{urec}) {
                        print "#####\n" if $verbose;
                        print $h{$oldpos}{$r}{mrec}."\n".$h{$oldpos}{$r}{urec}."\n" if $verbose;
                        print "#####\n" if $verbose;
                        print O $h{$oldpos}{$r}{mrec}."\n".$h{$oldpos}{$r}{urec}."\n";
                    }
                }
                # remove the items for the previous loci
                delete $h{$oldpos};
                print "deleted the hash item for locus $oldpos\n" if $verbose;
            }

            if ($a[1] =~ m/u/) { # unmapped read
                $h{$pos}{$a[0]}{urec} = $_;
                $oldpos = $pos;
                print "unmapped.. inserted into hash for $pos:$a[0] urec..\n" if $verbose;
            } elsif ($a[1] =~ m/U/ && $a[4]>0 && (m/XT:A:U/ || !m/XT:A/)) { # unique map
                $h{$pos}{$a[0]}{mrec} = $_;
                $oldpos = $pos;
                print "mapped.. inserted into hash for $pos:$a[0] mrec..\n" if $verbose;
            }
        }
    }
    # process the last loci
    for my $r (keys %{$h{$oldpos}}) {
        if (exists $h{$oldpos}{$r}{mrec} && exists $h{$oldpos}{$r}{urec}) {
            print O $h{$oldpos}{$r}{mrec}."\n".$h{$oldpos}{$r}{mrec}."\n";
        }
    }
    close F;
    close O;
}
