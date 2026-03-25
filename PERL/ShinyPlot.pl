#!/usr/bin/perl

# Run a sliding window analysis of haplotype divergence within a specified window size

# Usage: perl ShinyPlot.pl <path/to/Sliding_4R_outfile> <window_size> <step-size> <output-directory> <scaled? (yes/no)>

print "Usage: perl Slide_compare.pl <path/to/haplotypes-file> <window_size> <step-size> <output-directory>\n" if @ARGV < 5;

#use warnings;

#use strict;

# Declare global variables

my $Chr;

my $LineCount;

my @Positions;

my @DiffList;

my %PopulationHash;

my %HaplotypeHash;

my @RefHaplotypeList;

my $HaplotypeLen;

#my ($inFile, $windowSize, $stepSize, $outdir) = @ARGV;

#print "($inFile, $windowSize, $stepSize, $outdir)\n";

my ($inFile, $windowSize, $stepSize, $outdir, $scaled) = @ARGV if @ARGV == 5;

# open haplotypes file

open(DATA, $inFile) || die "Can't find infile $inFile\n";

while($L = <DATA>) {

  chomp($L);

  next if $L =~ /^snps/;

  if($L =~ /^sites/) {

    ($lineIdentifier, $sequence, $sites) = split(/\t/, $L);

    $SitesHash{$sequence} = $sites

  }

  else {

    ($strain, $pop, $sequence, $haplotypeString) = split(/\t/, $L);	# write haplotype data to a local array 

    $PopulationHash{$strain} = $pop;

    next if length($haplotypeString) < $ARGV[1] * 10;				# skip if number of variant sites is small (< ~50 sliding windows)

    $HaplotypeHash{$strain}{$sequence} = $haplotypeString;    				# write haplotype array to a hash keyed by strain ID and sequence 

  }

}

mkdir "$outdir" || die "Can't create output directory $!\n";

foreach my $TestStrain (sort {$a cmp $b} keys %HaplotypeHash) {

  foreach $sequence (sort {$a cmp $b} keys %{$HaplotypeHash{$TestStrain}}) {

    my @WindowStarts;

    my $sites = $SitesHash{$sequence};
    
    @Positions = split (/ /, $sites);

    my $NumSites = @Positions;

    print "$TestStrain\t$NumSites\n";

    open (OUT, '>', "$outdir/$sequence.$TestStrain.diffs") || die "Can't create outfile\n";

    for(my $j=0; $j<= $NumSites - $windowSize; $j += $stepSize) {

#      print "$Positions[$j]\n";

      push @WindowStarts, $Positions[$j]

    }

    my $NumWindows = @WindowStarts;

    print OUT "sites\t";

    print OUT join ("\t", @WindowStarts), " clade\n";		# print header line for R input

    my $TestHaplotype = $HaplotypeHash{$TestStrain}{$sequence};

    Compare_2_others($TestStrain, $sequence, $TestHaplotype);

    @WindowStarts = ();

    close OUT

  }

}

sub Compare_2_others {

  my $i;

  my($TestStrain, $sequence, $TestHaplotype) = @_;

  foreach my $comparatorStrain (sort {$a cmp $b} keys %HaplotypeHash) {

    next if $comparatorStrain eq $TestStrain;

    my $compHaplotype = $HaplotypeHash{$comparatorStrain}{$sequence};

    slidingWindows($comparatorStrain, $TestHaplotype, $compHaplotype);

  }

}

sub slidingWindows {

  my $j;

  my($comparatorStrain, $TestHaplotype, $Haplotype) = @_;

  my @WindowSNPsTotal = ();

  for($j=0; $j<= length($TestHaplotype) - $windowSize; $j += $stepSize) {

    my $TestHaplotypeWindow = substr($TestHaplotype, $j, $windowSize);

    my $HaplotypeWindow = substr($Haplotype, $j, $windowSize); 

    my $SNPsTotal = ( $TestHaplotypeWindow ^ $HaplotypeWindow ) =~ tr/\0//c;      

    my $actualWindowLength = $Positions[$j+$windowSize] - $Positions[$j+1];

    if($scaled eq 'yes') { 

      push @WindowSNPsTotal, $SNPsTotal/$actualWindowLength;

    }

    else {

      push @WindowSNPsTotal, $SNPsTotal/$windowSize;

    }

  }

  print OUT "$comparatorStrain\t";

  print OUT join ("\t", @WindowSNPsTotal), "\t$PopulationHash{$comparatorStrain}\n";

}       

close DATA
   

  

  
