#!/usr/bin/perl

##
# bwa_log_create.pl
#
# Description:
#       recreates a deleted or non-existent bwa.log file
#
# Version:
#       $Id: bwa_log_recreate.pl,v 1.4 2011/12/03 20:49:16 psm13 Exp psm13 $
#
# Revisions:
#       $Log: bwa_log_recreate.pl,v $
#       Revision 1.4  2011/12/03 20:49:16  psm13
#       set up tallying the different mapping types (uniq, non-uniq,
#       unmapped) and printing in the bwa.log format, I think
#
#       Revision 1.3  2011/12/03 20:44:00  psm13
#       gets the read groups from the header
#
#       Revision 1.2  2011/12/02 21:33:45  psm13
#       based on Lixing Yang's cntread_rg.pl
#
#       Revision 1.1  2011/12/02 21:32:09  psm13
#       Initial revision
#
# @Author: Psalm Mizuki (psm3426@gmail.com)
##

# cntread_rg.pl
# count reads per read group
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA
# Email: lixing_yang@hms.harvard.edu
# i.e. samtools view -X bamfile |perl scripts/cntread_rg.pl

use strict;

#argument check
if(scalar(@ARGV) != 1) {
    die "Usage: bwa_log_create.pl <BAM file or output prefix>";

}
my $input = $ARGV[0];
#make output file name, it's bamfile name (minus [.sorted].bam) + .bwa.log
my ($outfile, $bam);
if($input =~ /\.bam$/) {
	$outfile = substr($input, 0, -4);
	if($outfile =~ /\.sorted$/) {
		$outfile = $`; #gets the part right before the .sorted
	}
	$outfile .= ".bwa.log";
    $bam = $input;
} else {
	$outfile = $input . '.bwa.log';
	system "rm $outfile" if (-e $outfile);
    if(-e "$input.sorted.bam") {
        $bam = "$input.sorted.bam";
    } else {
        $bam = "$input.bam";
    }
}


my (%count, $total, $totUniq, $totNonUniq, $totUnmap);

#get the read groups from the BAM file header
open(PIPE, "samtools view -H $bam |") or die "Can't open pipe: $!\n";
my @pipe = <PIPE>;
close(PIPE);

foreach my $i (@pipe) {
    chomp($i); #get rid of the trailing newline
    if($i =~ /^\@RG/) {
        my @split = split(/\t/, $i); #isolate the ID:<read group>
        my $id = $split[1];
        @split = split(/ID:/, $id); #get rid of the ID: part
        $id = pop(@split);
        #print "id: $id\n";
        #assign to an array to hold total, uniq-map, non-uniq-map, and unmapped
        $count{$id} = [0, 0, 0, 0];
    }
}

if($input =~ /\.bam$/) {
	#open BAM file
	open(PIPE, "samtools view -X $bam |") or die "Can't open pipe: $!\n";
	while(my $line = <PIPE>) { #parse line by line
		chomp($line); # get rid of trailing newline
		my @data = split(/\t/, $line);
		for(my $i=11; $i<@data; $i++) { #skip the mandatory fields at beginning
			if ($data[$i] =~ m/RG:Z:/) { #get the read group of the read
				my $id = $';#grabs what came after the match (RG:Z:)
				#print "id: $id\n";
				++$count{$id}[0]; #increment the total for the read group
				++$total;
				if($line =~ /XT:A:/) { #if using XT tags
					if($line =~ /XT:A:U/) { #increment uniquely mapped
						++$count{$id}[1];
						++$totUniq;
					} elsif($line =~ /XT:A:R/) { #increment multiply mapped
						++$count{$id}[2];
						++$totNonUniq;
					} else { #anything else is unmapped
						++$count{$id}[3];
						++$totUnmap;
					}
				} else {
					if($data[1] =~ /u/) { #read is unmapped (bitwise flag)
						++$count{$id}[3];
						++$totUnmap;
					} elsif($data[4] > 0) { #unique reads have mapping qual > 0
						++$count{$id}[1];
						++$totUniq;
					} else { #multiply mapped reads
						++$count{$id}[2];
						++$totNonUniq;
					}
				}
				#print "\$count{$id} = " . $count{$id} . "\n";
			}
		}
	}
	close(PIPE);
} else {
	foreach my $id (keys %count) {
		my $detailfile = $input . '.rg' . $id . '.bwa_detail';
        #print "$detailfile\n";
		my $newline;
		open FILE, "<$detailfile";
		while($newline = <FILE>) {
			chomp $newline;
			my @data = split (/\t/, $newline);
			$count{$id} = \@data;
			$total = $data[0];
			$totUniq = $data[1];
			$totNonUniq = $data[2];
			$totUnmap = $data[3];
		}
		close FILE;
		#system "rm $detailfile" if (-e $detailfile);
	}
}

#print out everything
#print header and total counts and percentages
open(OUT, ">", "$outfile") or die "Can't open $outfile for writing: $!\n";
print OUT "\t#total\t#uniq\t#non-uniq\t#unmappable\t%uniq\t%non-uniq\t";
print OUT "%unmappable\ntotal:\t$total\t$totUniq\t$totNonUniq\t$totUnmap\t";
print OUT &percent($totUniq, $total) . "\t" . &percent($totNonUniq, $total);
print OUT "\t" . &percent($totUnmap, $total) . "\n";

#print out values for each read group
foreach my $key (keys %count) {
    my @tmp = @{$count{$key}}; #grab array from hash
    #print numerical values first
	print OUT "$key\t$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t";
    #print percentages
    print OUT &percent($tmp[1], $tmp[0]) . "\t" . &percent($tmp[2], $tmp[0]);
    print OUT "\t" . &percent($tmp[3], $tmp[0]) . "\n";
}

sub percent() {
    my $num = shift(@_);
    my $denom = shift(@_);

    if(!$denom) {
        return 0.00
    }
    return ($num/$denom * 100);
}
