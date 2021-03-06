#!/usr/bin/perl

# Run a sliding window analysis of haplotype divergence within a specified window size

# Reads outfiles generated by Create_slide4R.pl script

# Usage: perl Slide_compare.pl <path/to/Sliding_4R_outfile> <window_size> <step-size> <output-directory>

print "Usage: perl Slide_compare.pl <path/to/Sliding_4R_outfile> <window_size> <step-size> <output-directory>\n" if @ARGV != 4;

use warnings;

use strict;

# Declare global variables

my $Chr;

my $LineCount;

my @Positions;

my @DiffList;

my %PopulationHash;

my %HaplotypeHash;

my @RefHaplotypeList;

my $HaplotypeLen;

my ($inFile, $windowSize, $stepSize, $outdir) = @ARGV;



open(DATA, $inFile) || die "Can't find infile $inFile\n";

# Grab chromosome number from filename 

if($inFile =~ /(Chr\d)_.+/) {

  $Chr = $1;

#  print "$Chr\n";

}


while(my $L = <DATA>) {

  $LineCount ++;

  chomp($L);

  if($LineCount == 1) {	

    @Positions = split(/ /, $L); 		# write SNP positions to a global array

    unshift @Positions, "  ";

#    print "@Positions\n"

  }

  else {

    my @Haplotype = split(/ /, $L);		# write haplotype data to a local array 

    my $Strain = shift @Haplotype;		# grab Strain ID from front of array

    $PopulationHash{$Strain} = $Haplotype[-2];	# grab population ID from penultimate position in array

    splice(@Haplotype, -2);			# remove last two elements of array

    $HaplotypeLen = @Haplotype - 1;			# determine number of elements in haplotype 

    $HaplotypeHash{$Strain} = [@Haplotype];    	# write haplotype array to a hash keyed by strain ID 

  }

}

mkdir "$outdir" || die "Can't create output directory $!\n";

foreach my $TestStrain (keys %HaplotypeHash) {

  my @WindowStarts;

  my $NumSites = @Positions;

  print "$TestStrain\n";

  open (OUT, '>', "$outdir/$Chr.$TestStrain.diffs") || die "Can't create outfile\n";

  for(my $j=1; $j<= $NumSites - $windowSize; $j += $stepSize) {

    push @WindowStarts, $Positions[$j]

  }

  my $NumWindows = @WindowStarts;

  print OUT "\t";

  print OUT join ("\t", @WindowStarts), " clade\n";		# print header line for R input

  @RefHaplotypeList = @{$HaplotypeHash{$TestStrain}};

  Compare_2_Others($TestStrain);

  @WindowStarts = ();

  close OUT

}

sub Compare_2_Others {

  my $i;

  my($TestStrain) = @_;

  foreach my $Strain (keys %HaplotypeHash) {

    @DiffList = ();

    next if $Strain eq $TestStrain;

    my @HaplotypeList = @{$HaplotypeHash{$Strain}};

    for($i = 0; $i <= $HaplotypeLen; $i++) {

#      my($k, $l) = ($RefHaplotypeList[$i], $HaplotypeList[$i]);

#      print "$k, $l\n";

      push @DiffList, abs($RefHaplotypeList[$i] - $HaplotypeList[$i]);

    }

    slidingWindows($Strain);

  }

}

sub slidingWindows {

  my $j;

  my($Strain) = @_;

  my @WindowSNPsTotal = ();

  for($j=0; $j<= @DiffList - $windowSize; $j += $stepSize) {

    my $SNPsTotal = 0;

    foreach my $Value (@DiffList[$j..$j+$windowSize-1]) {

      $SNPsTotal += $Value

    }

    push @WindowSNPsTotal, $SNPsTotal/$windowSize;

  }

  print OUT "$Strain\t";

  print OUT join ("\t", @WindowSNPsTotal), "\t$PopulationHash{$Strain}\n";

}       

close DATA
   

  

  
