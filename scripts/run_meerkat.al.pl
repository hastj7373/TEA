#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

#get command line options then after they're gone, do an argument check
my %options;
getopts('astp:r:', \%options);

#args: BAM file, step (of the Meerkat program), output file prefix
if(scalar(@ARGV) < 2) {
    my $usage = "Usage: perl run_meerkat.al.pl [-a -p <matched BAM> -r <step> ";
    $usage .= "-s -t] <BAM> <reference genome>\nOptions:\n";
    $usage .= "\t-a: BAM was mapped using an aligner other than BWA\n";
    $usage .= "\t-p <matched BAM>: if current BAM is part of a tumor-normal ";
    $usage .= "pair, list\n\t\t\t  the location of the normal BAM\n";
    $usage .= "\t-r <step>: rerun the given step\n";
    $usage .= "\t-s: short reads (<50bp)\n\t-t: run only first two steps\n";
	die $usage;
}

#assign information from command line
my $short_reads = $options{s};
my $three_steps = $options{t};
my $other_aligner = $options{a};
my $step_start = 0;
my $step_end = 4;
if(defined($options{r})) {
    $step_start = $options{r};
}

my $samplebam;
my $ref;
my $th_num;

if(scalar(@ARGV) == 3) {
    ($samplebam, $ref, $th_num) = @ARGV;
} elsif(scalar(@ARGV) == 2) {
    ($samplebam, $ref) = @ARGV;
    $th_num = 1;
}

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

#programs and locations being used over and over
my $bwa = "$hash{'bwa_home'}";
if(($ref eq "human_g1k_v37")||($ref eq "hs37d5")) {
    $bwa = "/home/sl279/BiO/World/usr/local/bin/";
}
my $fasta_loc = "$hash{'fasta_location'}";
my $samtools = "$hash{'samtools_home'}";
my $blast = "$hash{'blast_home'}";
my $meerkat_home = "$hash{'meerkat_home'}";
my $tea_home = "$hash{'tea_home'}";

my $pre_process = "$meerkat_home/pre_process.pl";
my $meerkat = "$meerkat_home/meerkat.pl";
my $rmEmptyFqs = "$tea_home/scripts/rmEmptyFqs.pl";
my $bwaLogCreate = "$tea_home/scripts/bwa_log_create.pl";
my $generateBlacklist = "$tea_home/scripts/generateBlacklist.pl";

#derive other needed values#
my $bam_base_name = substr($samplebam, 0, -4); #file name minus .bam
if($bam_base_name =~ /\.sorted$/) { #take off .sorted too if it's there
    $bam_base_name = $`; #gets the part right before the .sorted
}
print "bam_base_name: $bam_base_name\n";

my $bwalogfile = "$bam_base_name.bwa.log";
my $blackfile = "$bam_base_name.blacklist.rg";
my $output = $bam_base_name; #output directory for fastqs

#create job name base, full name depends on program step
#isolate the BAM name
my $rind = rindex($bam_base_name, '/', (length($samplebam)-2)); 
my $jobname_base = substr($bam_base_name, $rind+1); #only the base name
print "job_name_base: $jobname_base\n";

#create bsub command and shell script for each step
for(my $step = $step_start; $step < $step_end; ++$step) {
    #create the jobname
    my $jobname = "$jobname_base-M$step";

    #the contents of the shell script file being run
    my $script = "#!/bin/bash\n\n";
    $script .= "echo \"=======================$step=======================\"";
    $script .= ">&2\n";

    if($step == 0) {
        #if there isn't a BWA log file, then create one
        unless(-e $bwalogfile) {
            $script .= "echo \"creating a BWA log file\">&2\n";
            $script .= "perl $bwaLogCreate $samplebam\n";
        }
        unless(-e $blackfile) {
            $script .= "echo \"creating a blacklist if necessary\">&2\n";
            $script .= "perl $generateBlacklist $bwalogfile";
        }
        
    } elsif($step == 1) {
	$script .= "perl $pre_process -s 20";
        $script .= " -k 1500 -t $th_num -P is -b $samplebam -R $blackfile -I ";
        $script .= "$fasta_loc ";
        $script .= "-A $fasta_loc.fai -W $bwa ";
        $script .= "-S $samtools";
        if($ref eq "mm9") {
            $script .= " -q 5";
        }

        #add in step where empty .fq.gz files are deleted
        $script .= "\necho \"removing empty FASTQ files\">&2";
        $script .= "\nperl $rmEmptyFqs $output";
    } elsif($step == 2) {
	$script .= "perl $pre_process -s 20 ";
        $script .= "-k 1500 -t $th_num -P cl -b $samplebam  -R $blackfile -I ";
        $script .= "$fasta_loc ";
        $script .= "-A $fasta_loc.fai -W $bwa ";
        $script .= "-S $samtools ";
        if($ref eq "mm9") {
            $script .= "-q 5 ";
        }
        if($short_reads) {
            $script .= "-l 0"; #for Alice's short reads
        }
    } elsif($step == 3) {
    	$script .= "perl $meerkat -s 20 ";
        $script .= "-t $th_num -p 3 -o 1 -P dc -b $samplebam -R $blackfile ";
    	$script .= "-B$blast -S $samtools ";
        $script .= "-F $fasta_loc -W $bwa";
        if($other_aligner) {
            $script .= " -u 1";
        }
        if($short_reads) {
            $script .= " -l 0";
        }
    }

    #write contents of $script to a file and make it executable
    open(FILE, ">", "$output.$step.sh") or die "Cannot open file: $!\n";
    print FILE "$script\n";
    close(FILE);

    chmod(0755, "$output.$step.sh");

    if($three_steps && ($step == 2)) {
        last;
    }
}
