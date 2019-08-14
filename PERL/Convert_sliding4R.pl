#!/usr/bin/perl

# Reads Strain ID file ("strain-idfile") & CP file and 

# creates matrix to be used as input to Slide_compare.pl

die "\nUsage: perl Convert_sliding_4R.pl <strain.idfile> <CP_file>\n\n" if @ARGV < 2;

use warnings;

use strict;


# declare global variables

my %IDhash;

my @IDs;


# Open chromopainter idfile and generate strain => population hash

open(IDFILE, $ARGV[0]);

while(my $L = <IDFILE>) {

  chomp($L);

  my($ID, $Pop, $Include) = split(/ /, $L);		

  push @IDs, $ID;

  $IDhash{$ID} = $Pop;	# Store info on each strain's population affiliation

}

close IDFILE;


# open chromopainter file

open(CPFILE, $ARGV[1]);

if(@ARGV == 3) {

  open(OUT, '>', $ARGV[2]) || die "Can't create outfile\n";

}

elsif($ARGV[1] =~ /(chr\d+)/) {

  my $Chr = $1;

  open(OUT, '>', "$Chr"."_haplotypes.txt") || die "Can't create outfile\n";

}

else {

  die "You either need to specify outfile name as an argument, or the haplotype filename must start with 'chr1', 'chr2', etc.\n"  

}


my $lineNum = 0;

my $hapNum = 0;	# my in orig

while(my $L = <CPFILE>) {

  chomp($L);

  $lineNum++;

  printPos($L) if $lineNum == 3;

  $hapNum = printHaps($L, $lineNum, $hapNum) if $lineNum >3;

  
}

close CPFILE;

close OUT;


## print the SNP positions

sub printPos {

  my @Positions = split(/ /, $_[0]);

  shift @Positions;

  print OUT join (" ", @Positions), " clade strain\n";

}


## print the haplotypes (space separated, with clade ID info)

sub printHaps {

  my($hapLine, $lineNum, $hapNum) = @_;

  my @Haps = split(//, $hapLine);

  unshift @Haps, $IDs[$hapNum];

  print OUT join (" ", @Haps), " $IDhash{$IDs[$hapNum]} $IDs[$hapNum]\n";

  $hapNum++;

  return($hapNum)

}



    
