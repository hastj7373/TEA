#!/usr/bin/perl

##
# generateBlacklist.pl
#
# Description:
#       generates blacklist file of read groups that have fewer than 30% or
#       their reads uniquely mapped
#
# Version:
#       $Id: $
#
# Revisions:
#       $Log: $
# @Author: Psalm Mizuki (psm3426@gmail.com)
##

use strict;
use warnings;

if(scalar(@ARGV) != 1) {
    die "Usage: generateBlacklist.pl <bwa.log file>\n";
}

my $file = $ARGV[0];
#take off the "bam" at the end
my $outfile = substr($file, 0, -8);
$outfile .= ".blacklist.rg";

my @blacklist;
open(FILE, "<", "$file") or die "Can't open $ARGV[0]: $!\n";
while(my $line = <FILE>) {
    chomp($line);
    unless(($line =~ /time/) || ($line =~ /total/)) {
        my @parts = split(/\t/, $line); #need columns 0, 1, and 2
        if($parts[1]) {
           if(($parts[2]/$parts[1]) < 0.3) {
               push(@blacklist, $parts[0]);
           }
        } else {
           push(@blacklist, $parts[0]);
        }
    }
}
close(FILE);

if(scalar(@blacklist) > 0) {
    open(FILE, ">", "$outfile") or die "Can't open $outfile to write: $!\n";
    foreach my $i (@blacklist) {
        print FILE "$i\n";
    }
    close(FILE);
    print STDERR "blacklist generated\n";
} else {
    print STDERR "no blacklist file to generate\n";
}
