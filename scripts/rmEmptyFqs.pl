#! /usr/bin/perl

##
# rmEmptyFqs.pl
#
# Description:
#       removes empty fastqs and isize distributions
#
# Version:
#       $Id: bam2fastq.pl,v 1.3 2011/11/01 15:06:32 psm13 Exp psm13 $
#
# Revisions:
#       $Log: bam2fastq.pl,v $
#       Revision 1.3  2011/11/01 15:06:32  psm13
#       hard coded in the locations of the programs to run to make it
#       extensible to removing empty fastq files
#
#       Revision 1.2  2011/10/19 20:43:35  psm13
#       grabbed code from dirTraverse.pl
#
#       Revision 1.1  2011/10/19 20:28:14  psm13
#       Initial Revision
#
# @Author: Psalm Mizuki (psm3426@gmail.com)
##

use strict;
use warnings;

if(scalar @ARGV != 1) {
    die "perl rmEmptyFqs.pl <fq.gz & .is dir>";
}
my $home_dir = $ARGV[0];
die "Enter a readable directory\n" unless opendir DIR, "$home_dir";

#takes care of ending slash problem
if(!($home_dir =~ /.+\/$/)) {
    $home_dir .= "/";
}

&dir($home_dir);

#psuedocode:
#if(looking at a file) {
#    print out file name and size
#}
#else { #directory, so recurse
#    get files
#    for( every file in there) {
#        call the original method
#    }
#}

sub dir() {
    my $file =  $_[0]; #$file can be a directory or a file
    #base case: if a file, then print and get 
    if(!(-d $file)) {
        chop($file);
        if(-f $file) {
            my $index = rindex($file, "/", (length($file)-2));
            my $name = substr($file, $index+1);
            if($file =~ /\.(fq|fastq)\.gz$/) {
                #delete empty fq.gz file (20 bytes in size)
                my @stats = stat($file);
                if($stats[7] < 21) {
                    if(unlink($file) == 0) {
                       print "Could not delete $file\n";
                    }
                    #delete corresponding .is if it's empty & hasn't been
                    my $prefix = substr($file, 0, -8);
                    my $isize_file = "$prefix.is";
                    if(-f $isize_file) {
                        @stats = stat($isize_file);
                        if($stats[7] < 1) {
                            if(unlink($isize_file) == 0) {
                                print "Could not delete $isize_file\n";
                            }
                        }
                    }
                }
            }
        }
    }
    elsif(-d $file) {
        #recursive case: if a directory
        local(*DIR);
        opendir DIR, $file;
    
        my $files;
        my @file_list;
        while(defined($files = readdir(DIR))) {
            push(@file_list, $files);
        }
        for(my $i=0; $i<=$#file_list; ++$i) {
            if(($file_list[$i] eq ".")||($file_list[$i] eq "..")) {
                next;
            }
            my $next = "$file"."$file_list[$i]/";
            #print "$next\n";
            &dir($next);
        }
    }
}
