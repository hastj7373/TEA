#!/usr/bin/perl
# Alice E. Lee: ejalice.lee@gmail.com
##########
# revision
#########
# contig : generate contigs for the clipped sequences and/or mates of RAMs using CAP3
#                        contig generation for     
# contig2 : revising the contig generation to process germline events 

use strict;
no strict "refs";
use Getopt::Std;
use File::Basename;
use IO::File;
use IO::Compress::Gzip qw(gzip $GzipError);
use List::Util qw(min max);

my $cmd = basename($0);
my $blast_path = "/home/el114/download/ncbi-blast-2.2.26+/bin";

my $usage = qq{

  Usage:   $cmd <command> [options]

  Command: aln      run bwa and generate a bam (using sampe)

           rdm      generate (uqmap and) rdm 
                    (coordnates of reads with unique mapping or discordant mates)
                    and generate bams containing discordant pairs
                    and write the total discordant pair count in the output file *.dstat

					 cbam			extract mappings with clipped reads and make a bam

					 rbam			extract mappings for ram 

           rid       run rid  # call rid.sample in rid.r through run.rid.r 
                     input: sample sample_desc_file(includes isize) etc.

           uqm       extract unique maps for a bam or bam folder
                     and write the total cnt in the output file *.mstat
           
           ram       generate a ram file per repeat type from ref and repeat mappings

           bamram    generate a ram file from the ref and repeat bam files  
            
           sram      separte rams per each repeat type 

           isize		 extract isize and write to *.isize	

           cisize	   convert isize infor from Lixing's to mine and extract read length information 

           bsort     sort bams

           findmap	 find mapping of input reads in a bam folder/file 

           bam2fq   generate fastq from a bam
     
           contig   generate contigs for the clipped sequences and/or mates of RAMs using CAP3

           contig2   generate contigs for the clipped sequences and/or mates of RAMs using CAP3 using the final call set 
           map_consensus  <infile> <family> <map_subfmaily: 1|0> 
};

my $log_dir = ".";
my @chrl = (); 
my %chrl = ();
&get_chrl(\@chrl, \%chrl, "human");

if (@ARGV < 1) { die $usage }

my $cmd = shift;
my @functions = qw(aln isize cisize cbam bsort uqm rdm ram bamram sram rid bam2fq contig contig2 map_consensus make_subfamilyseq);

if ($cmd eq "map_consensus") {
	my ($infile, $family, $map_subfamily) = @ARGV;
	&map_consensus($infile, $family, $map_subfamily, 1);
	exit(1);
} 

if ($cmd eq "make_subfamilyseq") {
	my ($family) = @ARGV;
	&make_subfamilyseq($family);
	exit(1);
}

if (grep {$_ eq $cmd} @functions) { &$cmd(); }
else { die $usage; }

sub make_subfamilyseq
{
  my ($family, $map_subfamily, $rdbs) = @_;
	my ($cmd, $db, $diag);
  if ($family eq "L1") { 
    $db = "/home/el114/repeat_analysis/data/consensus/M80343.1.fasta";
	} elsif ($family eq "Alu") {
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
	} elsif ($family eq "Alu") {
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluYa5.fasta"); 
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluYb8.fasta"); 
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluYb9.fasta"); 
		push(@$rdbs, "/home/el114/repeat_analysis/data/consensus/AluSc.fasta"); 
	} 
}

# blast clipped and ram mate contigs to consensus sequences
# and annotate insertion size, inversion status and diagnostic features for L1
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
	while (<F>) {
		chomp; my @a=split(/\t/);
		$insertions{$a[$idx{id}]}{line} = $_;
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
	for my $id (keys %insertions) {
		print O2 "$insertions{$id}{line}";

		# summarize divs and extract subfamilies with the minimum div		
		my (@divs, @divpcts, @est_subfamily);
		my $mindiv = 9999;
		for my $sf (@asf) {
			my $ssf = $sf; # short sf removing M80343.1 prefix
			if ($family eq "L1" && $sf ne "M80343.1") { # skip L1 ref name
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
	if (!$verbose) { `rm $file.*.map` }
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

sub contig
{
	my $usage = qq{
  Usage:   contig [options] <input file>
 
  Options: -a STR path to the contig assembler (cap3) 
           -d STR basedir to the cluster file with read names and clipped sequences and disc. bam files 
					 -t STR tempdir for running cap3 (/scratch/el114/)
           -o STR outdir (.)
           -c     clipped sequence contig only
           -r     ram mate contig only
           -v     verbose

	};

	print @ARGV;
  my %opts = ();
  getopts("ha:d:t:o:vrcv", \%opts);
  if (@ARGV < 1 || exists($opts{h})) { die $usage };
  my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $conly = 0; $conly = 1 if (defined $opts{c});
  my $ronly = 0; $ronly = 1 if (defined $opts{r});
  my $assembler = "cap3"; $assembler = $opts{a} if (defined $opts{a});
  my $dir = "/files/CBMI/parklab/alee/ra2"; $dir = $opts{d} if (defined $opts{d});
  my $odir = "."; $odir = $opts{o} if (defined $opts{o});
  my $tdir = "/scratch/el114"; $tdir = $opts{d} if (defined $opts{d});

  my $infile = shift(@ARGV);
	my $ph = basename($infile); 
	if ($ph =~ /cancer/) { $ph = "cancer" } elsif ($ph =~ /normal/) { $ph = "normal" } else { $ph=""}

	open(F, "< $infile") || die $usage;

	my %can = (); # fill a hash for candidate regions
	my @can = <F>; chomp(@can);
	close F;

	my $header = $can[0]; 
	for ( my $i=1; $i<@can; $i++) {
		my @a = split(/\t/, $can[$i]);
		#print "$a[0]:$a[1]:$a[2]:$a[5] => $a[15]\n";
	  #$can{$a[0]}{$a[1]}{"$a[2]:$a[5]"}{type} = $a[15];  #{"ov0890"}{"chr15:84128261-84128542"} = "M_L1" 
	  #$can{$a[0]}{$a[1]}{"$a[2]:$a[5]"}{type} = $a[13];  #{"ov0890"}{"chr15:84128261-84128542"} = "M_L1" 
	  $can{$a[0]}{$a[1]}{"$a[2]:$a[5]"}{type} = $a[10];  #{"ov0890"}{"chr15:84128261-84128542"} = "M_L1" 
	  $can{$a[0]}{$a[1]}{"$a[2]:$a[5]"}{rec} = $can[$i]; # fill the record 
	}

		for my $s (keys %can) {
			if ($verbose) { print "processing $s\n"; }
			if (!$ronly) { &clipped_contig($s,  $ph, \%can, $dir, $verbose); }
			if (!$conly) { &ram_contig($s, $ph, \%can, $dir, $verbose); }
		}
		&make_contig(\%can, $odir, $ronly, $conly, $assembler, $verbose);

		my $outprefix = $odir."/".basename($infile);
		&write_contig($outprefix, \%can, $header);
}

sub contig2 {
  my $usage = qq {
  Usage:   contig2 [options] <input file>
 
  Options: -a STR path to the contig assembler (cap3) 
           -d path to RUN DIRECTORY (originally: STR basedir to the cluster file with read names and clipped sequences and disc. bam files)
           -o STR outdir (.)
           -m STR consensus mapping [L1|Alu|SVA] (default: no)
           -s subfamily mapping (default: no)
					 -t STR tempdir for running cap3 (/scratch/el114/)
           -c     clipped sequence contig only
           -r     ram mate contig only
           -v     verbose

  };

  my %opts = ();
  getopts("ha:d:m:o:vrcvst:", \%opts);
  if (@ARGV < 1 || exists($opts{h})) { die $usage };
  my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $conly = 0; $conly = 1 if (defined $opts{c});
  my $ronly = 0; $ronly = 1 if (defined $opts{r});
  my $map = "no"; $map = $opts{m} if (defined $opts{m});
  my $map_subfamily = 0; $map_subfamily = 1 if (defined $opts{s});
  my $assembler = "cap3"; $assembler = $opts{a} if (defined $opts{a});
  my $dir = "/files/CBMI/parklab/alee/ra2"; $dir = $opts{d} if (defined $opts{d});
  my $odir = "."; $odir = $opts{o} if (defined $opts{o});
  my $tdir = "/scratch/el114/$map"; $tdir = $opts{t} if (defined $opts{t});
	if (! -d $tdir) { system("mkdir -p $tdir") }
	

  my $infile = shift(@ARGV);
	# XXX: Joe: this doesn't seem to be right.  If I do not add the lines
	# below (copied from contig), then the code does not append _cancer or
	# _normal to the samples.
	my $ph = ""; # for compatibility
	my $ph = basename($infile); 
	if ($ph =~ /cancer/) { $ph = "cancer" } elsif ($ph =~ /normal/) { $ph = "normal" } else { $ph=""}

  open(F, "< $infile") || die "can't open $infile";

  my %can = (); # fill a hash for candidate regions
  my @can = <F>; chomp(@can);
  close F;

  my $header = $can[0];
	my @a = split(/\t/, $header); my %idx = ();
	for (my $i=0; $i<@a; $i++) {
		$idx{$a[$i]} = $i;
	}
	#while (my ($key, $value) = each %idx) {
	#print "$key => $value\n";
	#}
	
  for (my $i=1; $i<@can; $i++) {
    my @a = split(/\t/, $can[$i]);
		if (exists $idx{s} && exists $idx{e}) {
    	$can{$a[$idx{sample}]}{$a[$idx{chr}]}{"$a[$idx{s}]:$a[$idx{e}]"}{type} = $a[$idx{"rep.repeat"}];
    	$can{$a[$idx{sample}]}{$a[$idx{chr}]}{"$a[$idx{s}]:$a[$idx{e}]"}{rec} = $can[$i];
		} else {
			# for compatibility for somatic candidates,
    	$can{$a[$idx{sample}]}{$a[$idx{chr}]}{"$a[$idx{pram_start}]:$a[$idx{nram_end}]"}{type} = $a[$idx{"rep.repeat"}];
    	$can{$a[$idx{sample}]}{$a[$idx{chr}]}{"$a[$idx{pram_start}]:$a[$idx{nram_end}]"}{rec} = $can[$i];
		}
  }

    for my $s (keys %can) {
      if ($verbose) { print "processing $s\n"; }
      if (!$ronly) { &clipped_contig($s,  $ph, \%can, $dir, $verbose); }
      if (!$conly) { &ram_contig($s, $ph, \%can, $dir, $verbose); }
    }
    &make_contig2(\%can, $odir, $tdir, $ronly, $conly, $assembler, $verbose);

    my $outprefix = $odir."/".basename($infile);
		#print $outprefix;
		&write_contig2($outprefix, \%can, $header);

		if ($map ne "no") { 
			&map_consensus($outprefix.".contig", $map, $map_subfamily, $verbose) 
		} 
}


sub ram_contig
{
	my ($s, $ph, $href, $dir, $verbose) = @_;

	# extract ram read names for candidate regions
	my (%prn, %nrn, %pcrn, %ncrn); # caution: one read can be mapped to multiple clusters!
	for my $chr (keys %{$href->{$s}}) {
		my $hhref = $href->{$s}->{$chr};
		# get the ram read names
		my $file;
		if ($ph eq "") { $file = "$dir/$s/cluster/$s.$chr.cluster.detail"; } # na18506/7/8
		else { $file = "$dir/$s\_$ph/cluster/$s\_$ph.$chr.cluster.detail"; }
		if (! -e $file) { print "$s $chr does not have the cluster detail file!\n"; next; }
		open(F, "< $file") || die "can't open $file";
		print "opened $file\n";
		while (<F>) {
			chomp; my @a = split(/\t/);
			if (exists $hhref->{"$a[1]:$a[2]"}) {
				my @pr = split(",", $a[30]);
				my @nr = split(",", $a[31]);
				my @ppos = split(",", $a[24]);
				my @npos = split(",", $a[25]);
				my $region = "$chr:$a[1]:$a[2]";
				for (my $i=0; $i<@pr; $i++) {
					if ($pr[$i] =~ m/mu1$/ || $pr[$i] =~ m/mu2$/ || $pr[$i] =~ m/sc$/) { 
						$pcrn{$pr[$i]}{chr} = $chr; 
						$pcrn{$pr[$i]}{pos} = $ppos[$i]; # we need to look the seq and qual for the reads mapped to the neg
						#$pcrn{$pr[$i]}{region} = $region;
						push(@{$pcrn{$pr[$i]}{region}}, $region);
					} else { 
						$prn{$pr[$i]}{chr} = $chr; 
						$prn{$pr[$i]}{pos} = $ppos[$i];
						#$prn{$pr[$i]}{region} = $region;
						push(@{$prn{$pr[$i]}{region}}, $region);
					}
				}
				for (my $i=0; $i<@nr; $i++) {
					if ($nr[$i] =~ m/mu1$/ || $nr[$i] =~ m/mu2$/ || $nr[$i] =~ m/sc$/) { 
						$ncrn{$nr[$i]}{chr} = $chr;
						$ncrn{$nr[$i]}{pos} = $npos[$i];
						#$ncrn{$nr[$i]}{region} = $region;
						push(@{$ncrn{$nr[$i]}{region}}, $region);
					} else { 
						$nrn{$nr[$i]}{chr} = $chr;
						$nrn{$nr[$i]}{pos} = $npos[$i];
						#$nrn{$nr[$i]}{region} = $region;
						push(@{$nrn{$nr[$i]}{region}}, $region);
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

	# read a disc or cl.disc.bam to get the read seq and qual
	my ($bamf, $clbamf);
	if ($ph eq "") {
		$bamf = "$dir/$s/bam/$s.disc.bam";
		$clbamf = "$dir/$s/bam/$s.cl.disc.bam";
	} else {
		$bamf = "$dir/$s\_$ph/bam/$s\_$ph.disc.bam";
		$clbamf = "$dir/$s\_$ph/bam/$s\_$ph.cl.disc.bam";
	}
	print "looking for ".scalar keys(%prn).":".scalar keys(%nrn)." ram names from disc.bam\n";
	print "looking for ".scalar keys(%pcrn).":".scalar keys(%ncrn)." ram names from cl.disc.bam\n";

	&get_read_seq($href->{$s}, $bamf, \%prn, \%nrn, $verbose);
	&get_read_seq($href->{$s}, $clbamf, \%pcrn, \%ncrn, $verbose);
}

sub get_read_seq
{
	my ($href, $bamf, $phr, $nhr, $verbose) = @_;
	my ($seq, $qual);

	print "start reading $bamf for assembling repeat read sequences\n";
	open(F, "samtools view -X $bamf |");
	my (@a, $n, $s, @b);
	while (<F>) {
		@a = split(/\t/); $n = $a[0];
		if (exists $phr->{$a[0]}) {
			if ($phr->{$a[0]}->{chr} ne $a[2] || abs($phr->{$a[0]}->{pos}) ne $a[3]) {
			# pos ram mate
			for my $r (@{$phr->{$a[0]}->{region}})	{
				#my @b = split(":", $phr->{$a[0]}->{region});
				#if ($verbose) { print "prseq: found seq for $phr->{$a[0]}->{region}, $a[0], $a[9]\n"; }
				#push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{prseq}}, $a[9]);
				#push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{prqual}}, $a[10]);
				my @b = split(":", $r);
				if (!($a[1] =~ m/r/))  {
    			($a[9] = reverse $a[9]) =~ tr/gatcGATC/ctagCTAG/; # if flat does not have 'r' reverse complement
    			$a[10] = reverse $a[10];
				}
				if ($verbose) { print "prseq: found seq for $r, $a[0], $a[9]\n"; }
				push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{prseq}}, $a[9]);
				push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{prqual}}, $a[10]);
			}
			}
		}
		if (exists $nhr->{$a[0]}) {
			if ($nhr->{$a[0]}->{chr} ne $a[2] || abs($nhr->{$a[0]}->{pos}) ne $a[3]) {
			for my $r (@{$nhr->{$a[0]}->{region}})	{
				#my @b = split(":", $nhr->{$a[0]}->{region});
				#if ($verbose) {print "nrseq: found seq for $nhr->{$a[0]}->{region}, $a[0], $a[9]\n"; }
				#push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{nrseq}}, $a[9]);
				#push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{nrqual}}, $a[10]);
				my @b = split(":", $r);
				if ($a[1] =~ m/r/)  {
    			($a[9] = reverse $a[9]) =~ tr/gatcGATC/ctagCTAG/;
    			$a[10] = reverse $a[10];
				}
				if ($verbose) {print "nrseq: found seq for $r, $a[0], $a[9]\n"; }
				push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{nrseq}}, $a[9]);
				push(@{$href->{$b[0]}->{"$b[1]:$b[2]"}->{nrqual}}, $a[10]);
			}
			}
		}
	}
	print "done reading $bamf for assembling repeat read sequences\n";
	close F;
}

sub clipped_contig
{
	my ($s, $ph, $href, $dir, $verbose) = @_;

	# get clipped sequences and read names for each candidate region
	my @chrs = keys %{$href->{$s}};
	for my $chr (@chrs) {
		my $clfile;
		if ($ph eq "") { $clfile = "$dir/$s/cluster/$s.$chr.clipped"; }
		else { $clfile = "$dir/$s\_$ph/cluster/$s\_$ph.$chr.clipped"; }
		if (! -e $clfile) { print "$s $chr does not have the clipped sequence file $clfile!\n"; next; }
		open(F, "< $clfile") || die "can't open $clfile";
		print "opened $clfile\n";
		while (<F>) {
			chomp; my @a = split(/\t/);
			if (exists $href->{$s}->{$chr}->{"$a[1]:$a[2]"}) {
				if ($verbose) {	print "$s:$chr:$a[1]:$a[2]:"; }
				if ($a[5]>0) {
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pseq}}, $a[10]);
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pqual}}, $a[11]);
					if ($verbose) {	print "pseq".join(",", @{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{pseq}})."\n"; }
				} else {
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nseq}}, $a[10]);
					push(@{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nqual}}, $a[11]);
					if ($verbose) {	print "nseq".join(",", @{$href->{$s}->{$chr}->{"$a[1]:$a[2]"}->{nseq}})."\n"; }
				}
			}
		}
	}
}

sub write_contig2
{
  my ($outprefix, $href, $header) = @_;
  my $outfile = $outprefix.".contig";
 
  open(F, "> $outfile") || die "Can't create $outfile";
  print "creating $outfile\n";
  print F "$header\torientation\tpolyA\tpolyT\tpclipped\tnclipped\tprammate\tnrammate\n";

  for my $s (keys %{$href}) {
    for my $chr (keys %{$href->{$s}}) {
      for my $site (keys %{$href->{$s}->{$chr}}) {
        my @pos = split(":", $site);
        my $hhref = $href->{$s}->{$chr}->{$site};
				&is_polyA2($hhref);

				my $orientation = "NA";
				if ($hhref->{polyA} eq "polyA" && $hhref->{polyT} eq "-") { $orientation = "+" } 
				if ($hhref->{polyT} eq "polyT" && $hhref->{polyA} eq "-") { $orientation = "-" } 
				
        #print F "$s\t$chr\t$pos[0]\t$pos[1]\t$hhref->{type}\t";
        print F "$hhref->{rec}\t$orientation\t";
        print F "$hhref->{polyA}\t$hhref->{polyT}\t";
        print F "$hhref->{pcontig}\t$hhref->{ncontig}\t$hhref->{prcontig}\t$hhref->{nrcontig}\n";
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
	my ($href, $odir, $tdir, $ronly, $conly, $assembler, $verbose) = @_;

	my (@seq, @qual, @pcontig, @ncontig, @psinglet, @nsinglet, $cmd);
	for my $s (keys %{$href}) {
		for my $chr (keys %{$href->{$s}}) {
			for my $site (keys %{$href->{$s}->{$chr}}) {
				my $hhref = $href->{$s}->{$chr}->{$site};
				my $fprefix = "$tdir/$s.$chr.$site";
				$hhref->{pcontig} = "";
				$hhref->{ncontig} = "";
				$hhref->{prcontig} = "";
				$hhref->{nrcontig} = "";

				if (!$ronly) {
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

				if (!$conly) {
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

        $cmd = "rm $fprefix.*"; #print $cmd."\n";
        if (!$verbose) { system($cmd); }

				if ($verbose) {
					if (!$ronly) {
					print "$s:$chr:$site:pseq:".@{$hhref->{pseq}}."\n";
					print "$s:$chr:$site:pqual:".@{$hhref->{pqual}}."\n";
					print "$s:$chr:$site:pcontig:".$hhref->{pcontig}."\n";

					print "$s:$chr:$site:nseq:".@{$hhref->{nseq}}."\n";
					print "$s:$chr:$site:nqual:".@{$hhref->{nqual}}."\n";
					print "$s:$chr:$site:ncontig:".$hhref->{ncontig}."\n";
					}
					if (!$conly) {
					print "$s:$chr:$site:prseq:".@{$hhref->{prseq}}."\n";
					print "$s:$chr:$site:prqual:".@{$hhref->{prqual}}."\n";
					print "$s:$chr:$site:prcontig:".$hhref->{prcontig}."\n";

					print "$s:$chr:$site:nrseq:".@{$hhref->{nrseq}}."\n";
					print "$s:$chr:$site:nrqual:".@{$hhref->{nrqual}}."\n";
					print "$s:$chr:$site:nrcontig:".$hhref->{nrcontig}."\n";
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
	
	open(F, "samtools view $fname |") or die "Can't open $fname: $!";
	my $prefix = basename($fname);
	$prefix =~ s/.bam//;
	
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
  if ($organism eq "human") {
    for (my $i=1; $i<=22; $i++) { push(@$chrl_ref, "$i"); $chrl_href->{$i} = 1 };
    push(@$chrl_ref, "X"); $chrl_href->{"X"} = 1;
    push(@$chrl_ref, "Y"); $chrl_href->{"Y"} = 1;
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
#	my $outf = basename($bamf); 
	my $outf2 = basename($bamf);
	my $outf3 = basename($bamf);
#	$outf =~ s/\.bam/.softclips.bam/;
	$outf2 =~ s/\.bam/.softclips.consd.raw.bam/;
#	$outf = "$outdir/$outf";  
	$outf2 = "$outdir/$outf2";  
	$outf3 =~ s/\.bam/.softclips.consd.cpos.bz2/;
	$outf3 = "$outdir/$outf3";  

#	open(O, "| samtools view -bS - > $outf") || die "Can't create $outf";
	open(O2, "| samtools view -bS - > $outf2") || die "Can't create $outf2";
	open(O3, "| bzip2 - > $outf3") || die "Can't create $outf3";

	# print the bam header
  open(F, "samtools view -H $bamf |");
	my @header = <F>;
#	print O @header;
	print O2 @header;
  close F;

	open(F, "samtools view $bamf|") || die "Can't open $bamf";
	my ($cigar, $qual, $qual1, $qual2, $s1, $s2, $selected, $selected2, $cpos);
	while(<F>) {
		chomp; my @a=split(/\t/);
    if ($a[1] & 0x0004) { next } # unmapped ; there must not be this case for the clipped reads, but just in case
		$cigar = $a[5]; $qual = $a[10];
		$s1 = $cigar;
		$s2 = $cigar;
		$selected = $selected2= 0;	
    if ($s1 =~ /^(\d+)S\S+/) {
        if ($1 >= 5) {
        	$qual1 = substr($qual, 0, $1);
        	if ($qual1 =~ m/([^#]+)/) {
						 if (length($1) >= 5) { 
							$selected = 1;	
							if (! ($a[1] & 0x0010)) { 
								$selected2 = 1; 
								$cpos = $a[3];
							}
						}
					}
				}
    }
    if ($s2=~ /\D(\d+)S$/) {
      if ($1 >= 5) {
 	    	$qual2 = substr($qual, -$1);
 	     	if ($qual2 =~ m/([^#]+)/) { 
					if (length($1) >= 5 ) { 
						$selected = 1;	
						if ($a[1] & 0x0010) { 
							$selected2 = 1; 
						  $cpos = &get_cpos($a[3], $cigar, $qual2, -1);
						} 
					}
				} 
			}
    }
#		if ($selected) { print O "$_\n"; 
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
	$prefix =~ s/.raw.bam$//;
	system("samtools sort $outf2 $prefix; samtools index $prefix.bam; rm $outf2");
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

	# representative (median) for all read groups
	my $all_is_mean = int(&median(&get_hash_values(\%is_mean)));
	my $all_is_median = int(&median(&get_hash_values(\%is_median)));
	my $all_is_sd = int(&median(&get_hash_values(\%is_sd)));
	my $all_rl = int(&median(&get_hash_values(\%rl)));	

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

sub isize 
{
	my $usage = qq{
  Usage:   isize [options] <bam file/dir>

  Options: -v     verboase 
           -d     input is folder 
           -q STR queue name (all_2h)
           -m INT maximum isize (1000)
           -s INT sampling (500,000)

	};

	my %opts = ();
	getopts("hvdq:m:s:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $queue = "all_2h"; $queue = $opts{q} if (defined $opts{q});
	my $mis = 1000; $mis = $opts{m} if (defined $opts{m});
  my $nsampling = 500000; $nsampling = $opts{s} if (defined $opts{s});

  my $in = shift(@ARGV);

	if (defined $opts{d}) {
		my @files = <$in/*.rg*.bam>;
		my @jobs = (); 
		for my $f (@files) {
      my $base = basename($f); 
			$base =~ s/.bam//;
			my $jid = `bjobs -w | awk '/is.$base/ { print \$1 }'`;
      if (-e "$in/$base.isize" || $jid ne "" ) {  next; }

     	my $cmd = "bsub -q $queue -J is.$base -o $log_dir/is.$base.out \"ra.pl isize";
      if ($verbose) { $cmd = $cmd." -v"; }
			$cmd = $cmd." -m $mis -s $nsampling $f\" | perl -ne 'm/<(\\d+)>/; print \$1'";
       print "$cmd \n" if ($verbose);
       $jid = `$cmd`;
			 push(@jobs, $jid);
		}
		print "\n".scalar(@jobs)." jobs submitted\n";

		my @completed = ();
		print "checking job completion\n" if ($verbose);
		while (1) {
      @completed = ();
      for my $f (@files) {
        my $base = basename($f); $base =~ s/.bam//;
        if (-e "$in/$base.isize") { push(@completed, $base); }
      }
      print scalar(@completed)." jobs completed\n" if ($verbose);
      if (scalar(@completed) < scalar(@files)) {
        sleep(20);
      } else {
        last;
      };
    }
    print "done isize extraction.\n" if ($verbose);
    print "merging into ".basename($in).".isize\n" if ($verbose);
		&merge_isize($in, 1); # name.isize remove each isize file
	} else {
		my ($is, $is_sd) = &get_isize($in, $mis, $nsampling, $verbose);
	}
}

sub merge_isize 
{
	my ($d, $rm) = @_;
	my @files = <$d/*.isize>;
	my $outf = "$d/".basename($d).".isize";
	open(O, "> $outf") or die "can't create isize summary table $outf";
	my (@is, @is_sd);
	for my $f (@files) {
		$f =~ m/\.(rg\S+)\.isize/;
		my $rg = $1;
		open(F, "< $f") or die "can't read isize table $f";
    my $l = <F>; 
		chomp($l);
		my @a = split(/\t/, $l);
		push(@is, $a[0]);
		push(@is_sd, $a[1]);
		print O "$rg\t$l\n";
	}
	print O "all\t".&median(\@is)."\t".&median(\@is_sd)."\n";
	close O;

	if ($rm) {
		`rm @files`;
	}
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
  Usage:   bamram [options] <ref bam> <repeat bam files>

  Options: -v     verboase
           -q STR queue name (all_2h)
           -o STR out dir for ram.bam/ram.bz2

  };

	my %opts = ();
	getopts("hvq:o:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $queue = "all_2h"; $queue = $opts{q} if (defined $opts{q});
  my $outdir = "."; $outdir = $opts{o} if (defined $opts{o});

	my $refbam = shift(@ARGV);
	open(REF, "samtools view $refbam | ") or die "no reference bam: $refbam: $!";

	# load repeat mappings into a hash
	my %h = (); # record end(1/2):READNAME\tREPAT_NAME1:fr1,REPEAT_NAME2:fr2, ...
	my $base; 

	while (@ARGV) {
		my $rabam = shift(@ARGV);
		$base = basename($rabam); $base =~ s/(.*)_[1|2].ra.bam/$1/;

		my $end = $rabam;
		#$end =~ s/.*disc.sorted_(1|2).ra.bam/$1/;
		$end =~ s/.*disc_(1|2).ra.bam/$1/;
		print "reading bamfile: $rabam, end: $end\n";
		my ($mapcnt, $dupmapcnt) = (0, 0);
		open(RA, "samtools view $rabam | ") or die "no reapeat assembly  bam: $rabam: $!";
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

	my $ramf = "$outdir/$base.ram.bz2"; 
	my $rbamf = "$outdir/$base.ram.bam"; 

	open(O, "|bzip2 - > $ramf") or die "can't create $ramf";
	open(O2, "| samtools view -bS - > $rbamf") or die "can't create $rbamf";

  # print the bam header
	print "starting reading $refbam\n";
  open(F, "samtools view -H $refbam |");
  my @header = <F>;
  print O2 @header;
  close F;

	# matching read ids from the referent bam and generate a ram file and a bam only with rams
	open(F, "samtools view $refbam|" )or die "can't open $refbam";
	my $cnt=0;
	print "starting the id matching to generate bam from refbam ...\n";
  while (<F>) {
		#print;
    chomp;
    my @a = split(/\t/);
		$a[2] =~ s/chr//;
    if (!exists($chrl{$a[2]}) || $a[1] & 0x4 || !(m/XT:A:U/)) { next; }
    if ($a[1] & 0x40) {# end:1 
      if (exists $h{"2:$a[0]"}) {
				#print "end1:identified mate map:$_\n";
				my $str = "$a[0]:\"".$h{"2:$a[0]"}."\"";
				$str = substr($str, 0, min(length($str), 250));
				#$_ =~ s/$a[0]/$a[0]:$h{"2:$a[0]"}/;
				$_ =~ s/$a[0]/$str/;
        print O2 "$_\n";
				#print "$cnt: printed $str to $rbamf\n";
  	    if ($a[1] & 0x10) { 
         print O "$a[0]\t$a[2]\t-$a[3]\t".$h{"2:$a[0]"}."\n";
					#print "$cnt: printed $a[0]\t$a[2]\t-$a[3]\t".$h{"2:$a[0]"}." to $ramf\n";
       } else {
         print O "$a[0]\t$a[2]\t$a[3]\t".$h{"2:$a[0]"}."\n";
					#print "$cnt: printed $a[0]\t$a[2]\t$a[3]\t".$h{"2:$a[0]"}." to $ramf\n";
       }
        $cnt++;
      }
    } elsif ($a[1] & 0x80) { # end:2
      if (exists $h{"1:$a[0]"}) {
				#print "end2:identified mate map:$_\n";
				my $str = "$a[0]:\"".$h{"1:$a[0]"}."\"";
				$str = substr($str, 0, min(length($str), 250));
				#$_ =~ s/$a[0]/$a[0]:$h{"1:$a[0]"}/;
				$_ =~ s/$a[0]/$str/;
        print O2 "$_\n";
				#print "$cnt: printed $str to $rbamf\n";
       if ($a[1] & 0x10) { 
         print O "$a[0]\t$a[2]\t-$a[3]\t".$h{"1:$a[0]"}."\n";
					#print "$cnt: printed $a[0]\t$a[2]\t-$a[3]\t".$h{"2:$a[0]"}." to $ramf\n";
       } else {
         print O "$a[0]\t$a[2]\t$a[3]\t".$h{"1:$a[0]"}."\n";
					#print "$cnt: printed $a[0]\t$a[2]\t$a[3]\t".$h{"2:$a[0]"}." to $ramf\n";
       }
        $cnt++;
      }
    }
  }
	close F;
	close O; 
	close O2;

  #system("echo $cnt > $outdir/$base.ram.cnt"); # will be map cnts
	print "done generating ram and rbam\n";
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

     	my $cmd = "bsub -R \"rusage[mem=3000]\" -q $queue -J ram.$base -o $log_dir/ram.$base.out \"ra.pl ram -v $f $in_repeat $out_dir \" | perl -ne 'm/<(\\d+)>/; print \$1'";
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

	#while (my ($key, $value) = each %rpr) {
	#	print "$key => $value\n";
	#}

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

	open(O2, "| samtools view -bS - > $rbam") || die "Can't create $rbam";

	# print the bam header
  open(F, "samtools view -H $bamf |");
	my @header = <F>;
	print O2 @header;
  close F;

	open(F, "samtools view $bamf|") || die "Can't open $bamf";
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
	open(O2, "| samtools view -bS - > $rbam") || die "Can't create $rbam";

	# print the bam header
  open(F, "samtools view -H $bamf |");
	my @header = <F>;
	print O2 @header;
  close F;

	open(F, "samtools view $bamf|") || die "Can't open $bamf";
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

sub uqm 
{
	my $usage = qq{
  Usage:   uqm [options] <bam file/folder> <out dir>

  Options: -v     verboase 
           -d     input is folder 
           -q STR queue name (all_2h)
	};
	my %opts = ();
	getopts("hvdq:", \%opts);
	if (@ARGV < 1 || exists($opts{h})) { die $usage }
	my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $queue = "all_2h"; $queue = $opts{q} if (defined $opts{q});

	my $in = shift(@ARGV);
	my $out_dir = shift(@ARGV);

  unless (-e $out_dir) { system("mkdir -p $out_dir") };
	my $jobcnt=0;

	if (defined $opts{d}) {
		my @files = <$in/*.rg*.bam>;
		for my $f (@files) {
      my $base = basename($f); 
			$base =~ s/.bam//;
      unless (-e "$out_dir/$base.mstat") {
      	my $cmd = "bsub -q all_2h -J uq.$base -o $log_dir/uq.$base.out \"ra.pl uqm -v $f $out_dir\"";
        print "$cmd\n" if ($verbose);
        system($cmd);
				$jobcnt++;
      }
		}
    #sleep(60*3);
		print $jobcnt." jobs submitted\n" if ($verbose);

		if ($jobcnt != 0) {
    my @completed = ();
    while (1) {
      print "checking job completion\n" if ($verbose);
      @completed = ();
      for my $f (@files) {
        my $base = basename($f); $base =~ s/.bam//;
        if (-e "$out_dir/$base.mstat") { push(@completed, $base); }
      }
      print scalar @completed." completed\n" if ($verbose);
      if (scalar(@completed) < scalar(@files)) {
        sleep(60);
      } else {
        last;
      };
    }
		}

		my @mfiles = <$out_dir/*.merged.gz>;
		if (@mfiles != 24) {
    	print "start merging\n" if ($verbose);
   	 &merge($out_dir, "uq");
		} else {
			print "merging already completed\n" if ($verbose);
		}
	} else {
		&uqm_bam($in, $out_dir, $verbose);
	}
}

sub uqm_bam
{
  my ($f, $out_dir, $verbose) = @_;
	my $cnt = 0;

  print "reading $f\n" if ($verbose);

  # create output uq per chromosome 
  my $base = basename($f); $base =~ s/.bam//;
  my %out_uq = ();
  my $fname_uq;
  for my $ch (@chrl) {
    $fname_uq = "$out_dir/$base.uq.$ch.gz";
    my $fh_uq = new IO::File "|gzip -c > $fname_uq";
    $out_uq{$ch} = $fh_uq;
  }

  open(F, "samtools view $f |");
  my @a=();
  while (<F>) {
    if (m/XT:A:U/) {
    	chomp; @a = split(/\t/);
      if ($a[1] & 0x0004 || !exists $out_uq{$a[2]}) { next; }
      if ($a[1] & 0x0010) {  # reverse strand mapping
      	print {$out_uq{$a[2]}} "-$a[3]\n";
      } else {
      	print {$out_uq{$a[2]}} "$a[3]\n";
			}
			$cnt++;
		}
  }
  close F;
	system("echo $cnt > $out_dir/$base.mstat"); # will be map cnts 
}

sub merge
{
 	my ($out_dir, $type) = @_;
	if ($type ne "uq" &&  $type ne "rdm") {	
		die("merge type uq or rdm are allowed");
	}
  my $id = basename($out_dir);
  for my $ch (@chrl) {
    my @files = <$out_dir/*.rg*.$type.$ch.gz>;
    my $out = "$out_dir/$id.$type.$ch.merged.gz";
    `zcat @files | gzip - > $out; rm @files`;
  }
}

sub rdm
{
  my $usage = qq{
  Usage:   rdm[options] <bam file/folder> <out dir>
 
  Options: -v     verboase 
           -d     input is folder 
           -b     generate a bam with discordant pairs
           -q STR queue name (all_2h)
  };

  my %opts = ();
  getopts("hvdcbq:", \%opts);
  if (@ARGV < 1 || exists($opts{h})) { die $usage }
  my $verbose = 0; $verbose = 1 if (defined $opts{v});
  my $bam = 0; $bam = 1 if (defined $opts{b});
  my $queue = "all_2h"; $queue = $opts{q} if (defined $opts{q});

  my $in = shift(@ARGV);
  my $out_dir = shift(@ARGV);

  unless (-e $out_dir) { system("mkdir -p $out_dir") };

  if (defined $opts{d}) {
    my @files = <$in/*.rg*.bam>;
    for my $f (@files) {
      my $base = basename($f); $base =~ s/.bam//;
      unless (-e "$out_dir/$base.dstat") {
	    	my $cmd = "bsub -q $queue -J rdm.$base -o $log_dir/rdm.$base.out \"ra.pl rdm";
				if ($verbose) { $cmd = $cmd." -v " }
        if ($bam) { $cmd = $cmd." -b " }
				$cmd = $cmd." $f $out_dir\"";
      	print "$cmd\n" if ($verbose);
				system($cmd);
      }
    }
    #sleep(60*3);
    my @completed = ();
    while (1) {
      print "checking job completion\n" if ($verbose);
      @completed = ();
      for my $f (@files) {
        my $base = basename($f); $base =~ s/.bam//;
        if (-e "$out_dir/$base.dstat") { push(@completed, $base); }
      }
      print "@completed completed\n" if ($verbose);
      if (scalar(@completed) < scalar(@files)) {
        sleep(60);
      } else {
        last;
      };
    }

    print "start merging\n" if ($verbose);
    &merge($out_dir, "rdm");

  } else {
  	&rdm_bam($in, $out_dir, $bam, $verbose);
  }
}

sub rdm_bam 
{
  my ($f, $out_dir, $bam, $verbose) = @_;
	my $cnt = 0;

  print "processing $f\n" if ($verbose);
 
	my ($is, $is_sd) = &get_isize($f, 1000, 500000, $verbose);
	print "isize: $is\t$is_sd\n" if ($verbose);

  # create rdm and discordant bams
  my $base = basename($f); $base =~ s/.bam//;

  my %out_rdm= ();

  my ($out_dbam, $fh_dbam, $f_rdm);

	if ($bam) {
  	$out_dbam = "$out_dir/$base.dbam";
 		$fh_dbam = new IO::File "|samtools view -bS - > $out_dbam";
	}

  for my $ch (@chrl) {
    $f_rdm= "$out_dir/$base.rdm.$ch.gz";
    my $fh_rdm= new IO::File "|gzip -c >> $f_rdm";
    $out_rdm{$ch} = $fh_rdm;
  }

	if ($bam) { # print the dbam header
  	open(F, "samtools view -H $f |");
    while (<F>) { print {$fh_dbam} $_; }
    close F;
	}

  open(F, "samtools view $f |");
  my @a=();
  while (<F>) {
   if (m/XT:A:U/) {
      chomp; @a = split(/\t/);
      if ($a[1] & 0x0004 || !exists $out_rdm{$a[2]}) { next; }

      if (&is_discordant($_, $is, $is_sd)) {
        if ($a[1] & 0x0010) {  # reverse strand mapping
          print {$out_rdm{$a[2]}} "-$a[3]\n";
        } else {
          print {$out_rdm{$a[2]}} "$a[3]\n";
        }
        if ($bam) {
					print {$fh_dbam} "$_\n";
				}
				$cnt++;
      }
    }
  }
  close F;
  system("echo $cnt > $out_dir/$base.dstat"); # will report the read counts
}

sub is_discordant 
{
  my ($map, $is, $is_sd) = @_;

  my @a=split(/\t/, $map);
  if ($a[1] & 0x8) { return (0) } # if the mate is unmapped, no discordant
  if ($a[6] ne "=") { return (1) }
  if (abs($a[8]) > $is + 3 * $is_sd || abs($a[8]) < $is - 3 * $is_sd) { return (1) }
  if ($a[3] < $a[7]) {
      if ($a[1] & 0x10) { return (1) }  # - and +
  } else {
      if (!($a[1] & 0x10)) { return (1) } # + and -
  }
  return (0);
}

sub load_isize
{
	open(F, "< @_[0]");
	<F>;
	chomp;
	my @a = split(/\t/);
	return($a[0], $a[1]);
}

sub get_isize
{
  my ($f, $max_isize, $nsampling, $verbose) = @_;
	my ($is, $is_sd);

	my $f_isize = basename($f);
	$f_isize =~ s/.bam/.isize/;
	$f_isize = dirname($f)."/".$f_isize;

	# if isize is already calculated, load them	
	if (-e $f_isize) { 
    print "reading isize from existing $f_isize\n" if ($verbose);
		($is, $is_sd) = &load_isize($f_isize);
    print "isize:$is\t$is_sd\n" if ($verbose);
		return($is, $is_sd);
	}

  open(F, "samtools view $f | ");
  my @a; my ($cnt, $sum, $sqsum) = (0, 0, 0);
  while (<F>) {
    chomp; @a = split(/\t/);
    if (($a[1] & 0x0002) && $a[8] != 0) {  # if bam flag says properl paired and isize is not 0
			$cnt++; $sum += abs($a[8]); 
			$sqsum += $a[8]*$a[8]; 
		}
    if ($cnt >= $nsampling) { 
			last; 
		}
  }
	close F;

 	$is = int($sum / $cnt); 
	$is_sd = int(sqrt($sqsum / $cnt - $is * $is));

	open (O, "> $f_isize");
	print O "$is\t$is_sd\n";
	close O;
  print basename($f).":$is:$is_sd\n" if ($verbose);
  print "writing isize to $f_isize\n" if ($verbose);
  print "isize:$is\t$is_sd\n" if ($verbose);

  return ($is, $is_sd);
}

sub aln
{
	my $usage = "$0 aln [-v] [-q queue:all_1d] <fastq folder>";
	if (@ARGV <1) { die $usage }
	
	my %opts = ();
	getopts("hvq:", \%opts);
	
	my $logdir = "/files/CBMI/parklab/alee/ra/gbm/log";
	#my $bamdir = "/files/CBMI/parklab/alee/ra/gbm/bam";
	my $bamdir = "/groups/pgp-secure/park-tcga/alee/gbm";
	my $tempdir = "/scratch/el114"; # space for sai and pp
	#my $ref = "/files/CBMI/parklab/alee/ra/data/assembly/hg18.allri.masked.fa";
	my $ref = "/files/CBMI/parklab/alee/ra/data/assembly/b36_w_chrM_woCHR.fa";
	
	my $l = 50;
	my $k = 2;
	my $n = 2;
	my $t = 50;
	
	my $queue = "all_1d"; $queue = $opts{q} if (defined($opts{q}));
	my $verbose = 0; $verbose = 1 if (defined($opts{v}));
	my $fqdir = shift(@ARGV);
	
	my $sample = basename($fqdir);
	my $new_sample = $sample;
	$new_sample =~ s/TCGA-//; $new_sample =~ s/^(.{2}-.{4}-.{3})(.*)/$1/;
	my @fq1 = <$fqdir/*_1.rg*.fastq.gz>; my $fq1_href = &rg_to_fq(\@fq1); 
	my @fq2 = <$fqdir/*_2.rg*.fastq.gz>; my $fq2_href = &rg_to_fq(\@fq2);
	
	my @rgl1 = keys %$fq1_href;
	my @rgl2 = keys %$fq2_href;
	my $rg_ref = &common(\@rgl1, \@rgl2);
	
	print "common rg: @$rg_ref\n" if ($verbose);
	
	my @submitted = ();

  unless (-e "$bamdir/$new_sample") { system("mkdir -p $bamdir/$new_sample") };
  unless (-e "$tempdir/$new_sample") { system("mkdir -p $tempdir/$new_sample") };
		
	for my $rg (@$rg_ref) {
		my $jid = `bjobs -w | awk '/ref.$new_sample.$rg / { print \$1}'`;

		# need to check the successful job completion from the log file
		if ($jid eq "" && -s "$bamdir/$new_sample/$sample.$rg.bam"  == 0 ) {
			print "bsub -R \"rusage[mem=4000]\" -o $logdir/ref.$new_sample.$rg.log -q $queue -J ref.$new_sample.$rg \"zcat $fq1_href->{$rg} | trim_fastq.pl $t | bwa aln -n $n -l $l -k $k $ref - > $tempdir/$new_sample/$sample.$rg.1.sai;
zcat $fq2_href->{$rg} | trim_fastq.pl $t | bwa aln -n $n -l $l -k $k $ref - > $tempdir/$new_sample/$sample.$rg.2.sai;rm $tempdir/$new_sample/pp.$sample.$rg.1 $tempdir/$new_sample/pp.$sample.$rg.2; 
mkfifo $tempdir/$new_sample/pp.$sample.$rg.1 $tempdir/$new_sample/pp.$sample.$rg.2; 
zcat $fq1_href->{$rg} | trim_fastq.pl $t > $tempdir/$new_sample/pp.$sample.$rg.1 & 
zcat $fq2_href->{$rg} | trim_fastq.pl $t > $tempdir/$new_sample/pp.$sample.$rg.2 & 
bwa sampe -s -o 1000 $ref $tempdir/$new_sample/$sample.$rg.1.sai $tempdir/$new_sample/$sample.$rg.2.sai $tempdir/$new_sample/pp.$sample.$rg.1 $tempdir/$new_sample/pp.$sample.$rg.2 | samtools view -bt b36_w_chrM_woCHR.fa.fai - > $bamdir/$new_sample/$sample.$rg.bam; 
rm $tempdir/$new_sample/pp.$sample.$rg.1 $tempdir/$new_sample/pp.$sample.$rg.2\"\n"; 

#   print "bsub -R \"rusage[mem=4000]\" -o $logdir/b36.$sample.$rg.log -q $queue -J ref.$sample.$rg \"zcat $fq1_href->{$rg} | trim_fastq.pl $t | bwa aln -n $n -l $l -k $k $ref - > $bamdir/$new_sample/$sample.$rg.1.sai; 
#    zcat $fq2_href->{$rg} | trim_fastq.pl $t | bwa aln -n $n -l $l -k $k $ref - > $bamdir/$new_sample/$sample.$rg.2.sai; 
#    rm pp.$sample.$rg.1 pp.$sample.$rg.2; 
#    mkfifo pp.$sample.$rg.1 pp.$sample.$rg.2; 
#   zcat $fq1_href->{$rg} | trim_fastq.pl $t > pp.$sample.$rg.1 & 
#   zcat $fq2_href->{$rg} | trim_fastq.pl $t > pp.$sample.$rg.2 & 
#   bwa sampe -s -o 1000 $ref $bamdir/$new_sample/$sample.$rg.1.sai $bamdir/$new_sample/$sample.$rg.2.sai pp.$sample.$rg.1 pp.$sample.$rg.2 | samtools view -bt b36_w_chrM_woCHR.fa.fai - > $bamdir/$new_sample/$sample.$rg.bam; 
#   rm pp.$sample.$rg.1 pp.$sample.$rg.2\"\n";
	
			push( @submitted, $rg); 
		}
	}
	
	print scalar(@submitted)." jobs submitted: @submitted\n";
}
	
sub rg_to_fq 
{
	my $ref = shift;
	my %h = ();
	for my $fq (@$ref) {
		$fq =~ m/\S+_[12].(rg\S+).fastq.gz/;
		$h{$1} = $fq;
	}
	return \%h;
}
	
sub common 
{
	my ($ref1, $ref2)  = @_;
	my @common = ();
	my %seen = ();
  foreach (@$ref1) { $seen{$_} = 1; }
  foreach (@$ref2) { if ($seen{$_}) { push(@common, $_); } }
	return \@common;
}
