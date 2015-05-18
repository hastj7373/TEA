#!/usr/bin/perl

use strict;
use Getopt::Std;
use Bio::DB::Sam;
use Bio::DB::Fasta;
use FindBin '$Bin';
#no warnings 'Bio::DB::Bam::Alignment';

my %opts = (k=>1, d=>3, p=>2, o=>0, q=>1, z=>1000000000, s=>20, m=>1, a=>1, l=>1, t=>1, R=>'null', P=>'all');
getopts("b:k:d:c:p:o:q:z:s:m:a:g:f:l:t:R:F:S:W:B:P:n:h", \%opts);

my $bamfile = $opts{b};
my $blacklist = $opts{k};
my $sd_cutoff_disc = $opts{d};
my $sd_cutoff_cl = $opts{c};
my $support_mps = $opts{p};
my $support_mpf = $opts{o};
my $support_reads = $opts{q};
my $sv_size_cutoff = $opts{z};
my $cut_sr = $opts{s};
my $remove_dup = $opts{m};
my $ad_align = $opts{a};
my $alt_map_max = $opts{g};
my $alt_map_max_clip = $opts{f};
my $clip = $opts{l};
my $threads_bwa = $opts{t};
my $black_rg = $opts{R};
my $reference_path = $opts{F};
my $samtools_path = $opts{S};
my $bwa_path = $opts{W};
my $blastall_path = $opts{B};
my $step = $opts{P};
my $line = $opts{n};
$sd_cutoff_cl = $sd_cutoff_disc unless (defined($opts{c}));
if (defined($opts{h}))
{
	&print_usage;
	die;
}
die "bam file not specified\n" unless (defined($bamfile));
if (defined($step))
{
	die "wrong step specified\nmust be one of the following: all|dc|cl|mpd|alg|srd|rf, default all\n" unless ($step =~ /all|dc|cl|mpd|alg|srd|ft|rf/);
}
#	control progress
my $discordant =1 if ($step eq 'all' or $step eq 'dc'); # extract discordant read pairs
my $cluster = 	1 if ($step eq 'all' or $step =~ 'cl'); # generate discordant read pair clusters
my $mpd = 	1 if ($step eq 'all' or $step eq 'mpd'); # call candidates from discordant clusters
my $alg = 	1 if ($step eq 'all' or $step eq 'alg'); # construct search space based on mpd candidates and align split reads to search space
my $srd = 	1 if ($step eq 'all' or $step eq 'srd'); # call SVs from pseudo split reads
my $filter =	1 if ($step eq 'all' or $step eq 'srd' or $step eq 'ft'); # filter outputs
my $blast_bp = 	1 if ($step eq 'all' or $step eq 'rf'); # blast break points reads to break points regions
my $cl =        0;
$cl =           1 if ($step eq 'cl1');
$cl =           2 if ($step eq 'cl2');
$cl =           3 if ($step eq 'cl3');

die "path to reference fasta files not specified\n" if ($alg and !(defined($opts{F})));

my $prefix = substr ($bamfile, 0, -4);
if ($prefix =~ /.sorted$/)
{
	$prefix = $`;
}
my $blacklistfile = $prefix.'.blacklist.gz' if ($opts{k});

my ($newline, %is);
# parse median and standard deviation of insert size from file
my $isinfofile = $prefix.'.isinfo';
die "$isinfofile file not exist\n" unless (-e $isinfofile);
open FILE, "<$isinfofile";
while ($newline = <FILE>)
{
	chomp $newline;
	if ($newline =~ /Read length/)
	{
		my ($trash, $rg) = split ("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'rl'} = $newline;
		if ($is{'rlu'})
		{
			$is{'rlu'} = $is{$rg}{'rl'} if ($is{$rg}{'rl'} > $is{'rlu'});
		}
		else
		{
			$is{'rlu'} = $is{$rg}{'rl'};
		}
		if ($is{'rld'})
		{
			$is{'rld'} = $is{$rg}{'rl'} if ($is{$rg}{'rl'} < $is{'rld'});
		}
		else
		{
			$is{'rld'} = $is{$rg}{'rl'};
		}
	}
	if ($newline =~ /Median/)
	{
		my ($trash, $rg) = split ("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'median'} = $newline;
	}
	if ($newline =~ /Standard deviation/)
	{
		my ($trash, $rg) = split ("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'sd'} = $newline;
		$is{$rg}{'isu'} = $is{$rg}{'median'} + $is{$rg}{'sd'}*$sd_cutoff_cl;
		$is{$rg}{'isu'} = int($is{$rg}{'isu'}) + 1;
		$is{$rg}{'isd'} = $is{$rg}{'median'} - $is{$rg}{'sd'}*$sd_cutoff_cl;
		$is{$rg}{'isd'} = int($is{$rg}{'isd'}) - 1;
		if ($is{'isu'})
		{
			$is{'isu'} = $is{$rg}{'isu'} if ($is{$rg}{'isu'} > $is{'isu'});
		}
		else
		{
			$is{'isu'} = $is{$rg}{'isu'};
		}
		#print "$rg\t$is{$rg}{'rl'}\t$is{$rg}{'median'}\t$is{$rg}{'sd'}\t$is{$rg}{'isu'}\t$is{$rg}{'isd'}\n";
	}
}

my $samtools_command;
if (defined($opts{S}))
{
	$samtools_path .= '/' unless ($samtools_path =~ /\/$/);
	$samtools_command = $samtools_path.'samtools';
}
else
{
	$samtools_command = 'samtools';
}
my $bwa_command;
if (defined($opts{W}))
{
	$bwa_path .= '/' unless ($bwa_path =~ /\/$/);
	$bwa_command = $bwa_path.'bwa';
}
else
{
	$bwa_command = 'bwa';
}
my $formatdb_command;
if (defined($opts{B}))
{
	$blastall_path .= '/' unless ($blastall_path =~ /\/$/);
	$formatdb_command = $blastall_path.'formatdb';
}
else
{
	$formatdb_command = 'formatdb';
}
my $blastall_command;
if (defined($opts{B}))
{
	$blastall_path .= '/' unless ($blastall_path =~ /\/$/);
	$blastall_command = $blastall_path.'blastall';
}
else
{
	$blastall_command = 'blastall';
}
# XXX: original code, didn't work for me
#$Bin =~ /scripts/;
#my $bin_path = $`.'bin/';
my $bin_path = $Bin . '/';
my $dre_command = $bin_path.'dre4';
my $scluster_command = $bin_path.'scluster2';

my $reference_db = Bio::DB::Fasta->new($reference_path) if ($alg or $blast_bp);

my $time0 = time;
my $local0 = localtime($time0);
print STDERR "$local0 Meerkat started\n";

&discord ($bamfile, $prefix, $blacklistfile, $clip, $dre_command, $sd_cutoff_disc, $isinfofile, $samtools_command, $remove_dup, $black_rg) if ($discordant);
&cluster ($prefix, $blacklistfile, \%is, $cut_sr, $ad_align, $alt_map_max, $alt_map_max_clip, $clip, $samtools_command, $cl, $scluster_command) if ($cluster);
&mpd (\%is, $prefix, $support_mps, $support_mpf) if ($mpd);
my $time = time;
my $local = localtime($time);
print STDERR "$local called events by read pairs\n" if ($mpd);
&alg($line, $prefix, $bwa_command, $reference_db, \%is, $cut_sr, $threads_bwa) if ($alg);
my $time = time;
my $local = localtime($time);
print STDERR "$local aligned split reads\n" if ($alg);
&srd($line, $prefix, \%is, $cut_sr, $sv_size_cutoff, $support_reads) if ($srd);
my $time = time;
my $local = localtime($time);
print STDERR "$local confirmed events by split reads\n" if ($srd);
&filter($prefix) if ($filter);
my $time = time;
my $local = localtime($time);
print STDERR "$local filtered events\n" if ($srd);
&blast_bp($line, $prefix, $reference_db, \%is, $cut_sr, $formatdb_command, $blastall_command) if ($blast_bp);
my $time = time;
my $local = localtime($time);
print STDERR "$local refined break points by local alignments\n" if ($blast_bp);

my $time1 = time;
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "Time used: $difference:$minutes:$seconds\n";

sub discord
{
	my $bamfile = shift;
	my $prefix = shift;
	my $blacklistfile = shift;
	my $clip = shift;
	my $dre_command = shift;
	my $sd_cutoff_disc = shift;
	my $isinfofile = shift;
	my $samtools_command = shift;
	my $remove_dup = shift;
	my $black_rg = shift;
	my $drelog = $prefix.'.dre.log';
	my $cl_bam = $prefix.'.cl.sorted.bam';
	my $disc_bam = $prefix.'.disc.bam';
	my $dup_file = $prefix.'.dup.bam';
	my $disc_cl_bam = $prefix.'.cl.disc.bam';
	my $dup_cl_file = $prefix.'.cl.dup.bam';
	
	die "bam file doesn't exist\n" unless (-e $bamfile);
	
print "going to run dre\n";
	if ($blacklistfile)
	{
		if ($remove_dup)
		{
			system "$dre_command -f -v -P -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file $bamfile >$drelog";
			system "$dre_command -f -v -P -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file $cl_bam >>$drelog" if ($clip);
		}
		else
		{
			system "$dre_command -f -v -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_bam $bamfile >$drelog";
			system "$dre_command -f -v -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_cl_bam $cl_bam >>$drelog" if ($clip);
		}
	}
	else
	{
		if ($remove_dup)
		{
			system "$dre_command -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file $bamfile >$drelog";
			system "$dre_command -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file $cl_bam >>$drelog" if ($clip);
		}
		else
		{
			system "$dre_command -f -v -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam $bamfile >$drelog";
			system "$dre_command -f -v -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam $cl_bam >>$drelog" if ($clip);
		}
	}
print "dre done\n";
	my $disc_sort = $prefix.'.disc.sorted';
	my $disc_cl_sort = $prefix.'.cl.disc.sorted';
	my $time = time;
	my $local = localtime($time);
	print STDERR "$local extracted discordant read pairs\n";
	
	system "$samtools_command sort -n $disc_bam $disc_sort";
	system "$samtools_command sort -n $disc_cl_bam $disc_cl_sort" if ($clip);
	system "rm $disc_bam";
	system "rm $disc_cl_bam" if ($clip);
	my $time = time;
	my $local = localtime($time);
	print STDERR "$local sorted discordant bam by name\n";
}
