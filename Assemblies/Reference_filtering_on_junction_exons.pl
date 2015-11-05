#!/usr/bin/perl 

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

#author: Shabhonam.caim@tgac.ac.uk####
my $infile1,
my $infile2,
my $infile3,
my $outfile1,
my $memory,
my $help;
my $queue;
my $flag=0;
my $jobname;
&GetOptions (
'1=s' => \$infile1,
'2=s' => \$infile2,
'3=s' => \$infile3,
'o=s' => \$outfile1,
'j=s' => \$jobname,
'q=s' => \$queue,
'm=s' => \$memory,
'h' => \$help,
) ;

#\t-j LSF job name
if ($help or (not $infile1) or not ($infile2) or not ($infile3)) {
print <<"TEXT";

USAGE:
-1 Ref gff 
-2 Junctions gff
-3 Sorted BAM file (for exon coverage) 
-o Output prefix name
-j Jobname
-m Memory (In MB eg:10000=10G (approx))
-q LSF queue (Test128,Test256,Prod128,Prod256)
-h help
TEXT

exit;
}


system("bsub  -J $jobname\_grepexon -q $queue  -o out_grepexon -R \"span[ptile=2]rusage[mem=2000] \" -n 2 \" grep -w \"exon\" $infile1 > exons.gff \"");

system("bsub -w 'ended(\"$jobname\_grepexon\")' -J $jobname\_coveragebedexon -q $queue  -o out_coverageexon -R \"span[ptile=4]rusage[mem=60000] \" -n 2 \"source bedtools-2.24.0; coverageBed -a exons.gff -b $infile3 > coverage_exon.gff\"");

system ("bsub -o out_addintrons -q $queue -J $jobname\_addintrons -R \"span[ptile=2]rusage[mem=5000]\" -n 2 \"source genometools-1.5.4;gt gff3 -addintrons  -retainids  $infile1  > Reference.added_introns.gff \"");

system("bsub -w 'ended(\"$jobname\_addintrons\")' -J $jobname\_grep -q $queue  -o out_grepintron -R \"span[ptile=2]rusage[mem=2000] \" -n 2 \" grep -w \"intron\" Reference.added_introns.gff > introns.gff \"");

system("bsub -w 'ended(\"$jobname\_grep\")' -J $jobname\_coveragebedintron -q $queue  -o out_coverageintron -R \"span[ptile=4]rusage[mem=$memory] \" -n 2 \"source bedtools-2.24.0; coverageBed -a introns.gff -b $infile2 > coverage_intron.gff\""); 

system("bsub -w 'ended(\"$jobname\_coveragebedexon\")' -J stats -q $queue -o out_stats \"bash bash_stats.sh \"  ");




 
 













