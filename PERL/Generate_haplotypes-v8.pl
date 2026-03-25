#!/usr/bin/perl

##############################
#
# Generate_haplotypes.pl
#
# written by Mark L. Farman
#
# Purpose: Read SNPcaller outfiles & BLAST reports, make a comprehensive lookup list of haplotypes
#
# Note: Must be re-run every time a new strain is added (unavoidable)
#
##############################

#use strict;

#use warnings;

die "Usage: perl Generate_haplotypes-v7.pl <strainList> <SNPFILE_DIR> <BLASTFILE-DIR> <HAPLOTYPES_OUTFILE> <REF> <CHR|CONTIG>\n" if @ARGV < 5;

($strainList, $snpDir, $blastDir, $haplOutfile, $refGenome, $seqPrefix) = @ARGV;

$haplOutfile = $haplOutfile.".txt" if $haplOutfile !~ /txt$/;

use FetchGenomeOrig;
use IsItRepeated;

# BLAST result files must be named according to the format: Query.Subject.BLAST
# The BLAST format used must be: -outfmt '6 qseqid sseqid qstart qend sstart send btop' 

### declare global variables

my %NumSNPs;

my %QueryCounts;

my %SubjectCounts;

my %TotalLen;


### MAIN

& READ_REF_GENOME;

& STRAIN_LIST;

& GRAB_BIALLELIC_VARIANTS;
  
open(HAPLOS, '>', $haplOutfile) || die "Can't open outfile\n";
 
& HEADER_REF_LINES;
       
& STRAIN_VARIANTS;

close HAPLOS;

($nmhaplosOutfile = $haplOutfile) =~ s/\.txt/\.complete\.txt/;
  
open(NMHAPLOS, '>', $nmhaplosOutfile) || die "Can't open outfile\n";

& PRINT_NO_MISSING_HAPLOTYPES;

close NMHAPLOS;

sub READ_REF_GENOME {

  print "Reading genome into a hash...\n";

  ### read in lengths and number of contigs 
  $refGenomeRef = FetchGenomeOrig::getSeqs($refGenome);
  %refGenome= %$refGenomeRef;

  foreach my $key (keys %refGenome) {
    ($newkey = $key) =~ s/.+?(\d+)$/$1/;
    $refGenome{$newkey} = $refGenome{$key};
    $len = length($refGenome{$newkey});
    $TotalLen{$newkey} = 'O' x $len;
    delete $refGenome{$key}
  }
}

sub STRAIN_LIST {
  print "Reading strains list...\n";
  open(SL, $strainList) || die "Can't open strain list\n";

  while($SL = <SL>) {

    chomp($SL);
    my ($strain, $pop, $incl) = split(/\t| +/, $SL);
    $StrainHash{$strain} = $pop

  }

  close SL,

}


# Loop through all SNP files and record all biallelic SNPs relative to Reference genome

sub GRAB_BIALLELIC_VARIANTS {

  $snpDir =~ s/\/$//;
  opendir(OUTFILES, $snpDir) || die "Can't open outfiles directory: $snpDir\n";
  @FilesList = readdir(OUTFILES);
  print "Identifying all relevant SNP loci\n";

  foreach $snpFile (@FilesList) {
    next if $snpFile !~ /out$/;
    ($refid, $test, $tail) = split(/_v_|_out/, $snpFile);
    next unless (exists($StrainHash{$test}));

    open(SF, "$snpDir/$snpFile") || die "Problem\n";
    print "\nReading SNPs in file: $snpFile\n";

    while($L = <SF>) {
      chomp($L);
      @Data = split(/\t/, $L);
      if(@Data == 7) {				# skip lines with repeat field
        ($Ref, $Alt, $RefPos, $AltPos, $RefNucl, $AltNucl, $Dir) = @Data;
        $Ref =~ s/.+?(\d+)$/$1/;
        next if ($RefNucl eq '-' || $AltNucl eq '-');
        $snp = $RefNucl.$AltNucl;
        next if $allAlleleHash{$Ref}{$RefPos} eq 'X';
        if($allAlleleHash{$Ref}{$RefPos} =~ /[AGTC]/) {
          $allAlleleHash{$Ref}{$RefPos} = 'X' if $allAlleleHash{$Ref}{$RefPos} ne $AltNucl
        }

        else {
          $allAlleleHash{$Ref}{$RefPos} = $AltNucl;
          $allSNPsHash{$Ref}{$RefPos} = $snp
        }
      }
      elsif(@Data == 6) {                          # run on files lacking strand information and skip lines with repeat field
        ($Ref, $Alt, $RefPos, $AltPos, $RefNucl, $AltNucl) = @Data;
        $Ref =~ s/.+?(\d+)$/$1/;
        next if ($RefNucl eq '-' || $AltNucl eq '-');
        $snp = $RefNucl.$AltNucl;
        next if $allAlleleHash{$Ref}{$RefPos} eq 'X';
    
        if($allAlleleHash{$Ref}{$RefPos} =~ /[AGTC]/) {
          $allAlleleHash{$Ref}{$RefPos} = 'X' if $allAlleleHash{$Ref}{$RefPos} ne $AltNucl
        }
 
        else {
          $allAlleleHash{$Ref}{$RefPos} = $AltNucl;
          $allSNPsHash{$Ref}{$RefPos} = $snp
        }
      }
    } 
    close SF
  }

  closedir OUTFILES;

# Create list of SNP positions on reference chromosomes

  $regex = qr/A|G|T|C/;

  open(SNPsList, '>', $haplOutfile.".snps");

  foreach my $contig (sort {$a cmp $b} keys %allAlleleHash) {

    @positions = sort {$a <=> $b} keys %{$allAlleleHash{$contig}};

    $numSNPs += @positions;

    foreach my $pos (@positions) {

      print SNPsList "$contig\t$pos\n";

      if($allAlleleHash{$contig}{$pos} =~ /($regex)/g) {

        push @{$biAlleleHash{$contig}}, $pos;
        push @{$biSNPsHash{$contig}}, $allSNPsHash{$contig}{$pos}  

      }
    }
  }

  close SNPsList;

  print "Number of SNPs: $numSNPs\n";

}


# write header lines (***currently prints out multiallelic sites***)

sub HEADER_REF_LINES {

  foreach my $seq (sort {$a cmp $b} keys %biAlleleHash) {
    my @sites = @{$biAlleleHash{$seq}};
    my @snps = @{$biSNPsHash{$seq}}; #    print "$seq\t$value\t@{$biAlleleHash{$seq}}\n";
    # edit next line for minimum length
    if(@sites >= 1) {
      print HAPLOS "sites\tsequence$seq\t";
      print HAPLOS join(" ", @sites), "\n";
      print HAPLOS "snps\tsequence$seq\t";
      print HAPLOS join(" ", @snps), "\n";
    }
  }
  foreach my $seq (sort {$a cmp $b} keys %biAlleleHash) {
    my @sites = @{$biAlleleHash{$seq}};
    $refHapl = '1' x @sites;
    print HAPLOS join("\t", ($refid, "ref", "sequence".$seq, $refHapl)), "\n" if @sites > 500;
  }
}


# Reopen SNP outfiles and call variants strain-by-strain

sub STRAIN_VARIANTS {
        
  $snpDir =~ s/\/$//;
  opendir(OUTFILES, $snpDir) || die "Can't open outfiles directory: $snpDir\n";
  @FilesList = readdir(OUTFILES);
  
  foreach $snpFile (@FilesList) {
    next if $snpFile !~ /out$/;
    ($ref, $test, $tail) = split(/_v_|_out/, $snpFile);
    next unless (exists($StrainHash{$test}));
    
    open(SF, "$snpDir/$snpFile") || die "Problem\n";
    print "Making SNP calls from file: $snpFile\n";
    $snpFileCount++;

    # reset SNPs hash for each strain;

    while($L = <SF>) {
      chomp($L);
      @Data = split(/\t/, $L);
      
      if(@Data == 7) {                          # skip lines with repeat field
        ($Ref, $Other, $RefPos, $OtherPos, $RefNucl, $OtherNucl, $Dir) = @Data;
        $Ref =~ s/.+?(\d+)$/$1/;
        next if ($RefNucl eq '-' || $OtherNucl eq '-');
        $SNPsHash{$test}{$Ref}{$RefPos} = $OtherNucl;
      }
      elsif(@Data == 6) {                          # skip lines with repeat field
        ($Ref, $Other, $RefPos, $OtherPos, $RefNucl, $OtherNucl) = @Data;
        $Ref =~ s/.+?(\d+)$/$1/;
        next if ($RefNucl eq '-' || $OtherNucl eq '-'); 
        $SNPsHash{$test}{$Ref}{$RefPos} = $OtherNucl;
      }

    }  
    close SF;

    $blastFile = $blastDir."/$ref.$test.BLAST";
    my $alignStringRef;
    $alignStringRef = IsItRepeated::REF($blastFile);
    CALL_VARIANTS($test, \%SNPsHash, $alignStringRef)


  }
  closedir OUTFILES;

}


sub CALL_VARIANTS {

  %Haplotypes = ();

  $count = 0;
  print "Calling variants...\n";
  ($test, $SNPsHashRef, $alignStringRef) = @_;
  %SNPsHash = %$SNPsHashRef;
  %alignString = %$alignStringRef;

  foreach $sequence (sort {$a <=> $b} keys %biAlleleHash) {
    next unless (exists($SNPsHash{$test}));
    foreach $pos (@{$biAlleleHash{$sequence}}) { 
      if(exists($SNPsHash{$test}{$sequence}{$pos})) {
        $Haplotypes{$sequence} .= 0;
        $haplLength = length($Haplotypes{$sequence});
        $Invariant{$sequence}{$haplLength} = 0 unless exists ($Invariant{$sequence}{$haplLength});
        $Invariant{$sequence}{$haplLength}++
      }
     
      elsif(substr($alignString{$sequence}, $pos-1, 1) == 1)  {
        $Haplotypes{$sequence} .= '1';
        $haplLength = length($Haplotypes{$sequence});
        $Invariant{$sequence}{$haplLength}++
      }

      else {
        $Haplotypes{$sequence} .= '9';
        $haplLength = length($Haplotypes{$sequence});
        $MissingHash{$sequence}{$haplLength} = 1
      }
      print "$pos/n" if $haplLength == 0
    }
  }

  @haplokeys = keys %Haplotypes;
  PRINT_HAPLOTYPES($test, \%Haplotypes);

}

sub CULL_INVARIANTS {

  print "#SNP files = $snpFileCount\n";

  foreach $seqid (sort {$a <=> $b} keys %Invariant) {
    foreach $site (sort {$a <=> $b} keys %{$Invariant{$seqid}}) {
      next if $site == 0;
      if($Invariant{$seqid}{$site} < 1) {
#        print "yes; $Invariant{$seqid}{$site}\n";
#        $MissingHash{$seqid}{$site} = 1;
      }
      elsif($Invariant{$seqid}{$site} == $snpFileCount) {
#        print "yes; $seqid, $site: $Invariant{$seqid}{$site}\n";
        $MissingHash{$seqid}{$site} = 1;
      }
    }  
  }
}


sub PRINT_HAPLOTYPES {

  print "Printing haplotypes\n\n";
  ($test, $HaplotypesRef) = @_;
  %Haplotypes = %$HaplotypesRef;
  @hapkeys = keys %Haplotypes;
  $pop = $StrainHash{$test};

  foreach $sequence (sort {$a <=> $b} keys %Haplotypes) {
    print HAPLOS join ("\t", ($test, $pop, $seqPrefix.$sequence, $Haplotypes{$sequence})), "\n" if length($Haplotypes{$sequence}) >= 500;
  }
  %Haplotypes = ()
}

sub PRINT_NO_MISSING_HAPLOTYPES {

  open(HAPLOS, $haplOutfile);

  while($H = <HAPLOS>) {
    chomp($H);

    if ($H =~ /^sites/) {
      my($prefix, $sequence, $sites) = split(/\t/, $H);
      $sequence =~ s/sequence//;
      @sites = split(/ /, $sites);
      foreach my $pos (sort {$b <=> $a} keys %{$MissingHash{$sequence}}) {
        splice(@sites, $pos-1, 1);
      }
      # edit next line depending on purpose
      if (@sites >= 1) {
        print NMHAPLOS "sites\t$seqPrefix$sequence\t";
        print NMHAPLOS join (" ", @sites), "\n";
      }
    }  

    elsif ($H =~ /^snps/) {
      my($prefix, $sequence, $snps) = split(/\t/, $H);
      $sequence =~ s/sequence//;
      @snps = split(/ /, $snps);
      foreach my $pos (sort {$b <=> $a} keys %{$MissingHash{$sequence}}) {
        splice(@snps, $pos-1, 1);
      }
      if (@sites >= 500) {
        print NMHAPLOS "snps\t$seqPrefix$sequence\t";
        print NMHAPLOS join (" ", @snps), "\n";
      }
    }  

    else {
      my($test, $pop, $sequence, $haplotype) = split(/\t/, $H);
      $sequence =~ s/sequence//;

      foreach my $pos (sort {$b <=> $a} keys %{$MissingHash{$sequence}}) {
        substr($haplotype, $pos-1, 1, '');
      }
      print NMHAPLOS join ("\t", ($test, $pop, $seqPrefix.$sequence, $haplotype)), "\n" if length($haplotype) >= 500;
    }
  }
  
  close HAPLOS;
  close NMHAPLOS;
}


sub ALREADY_COMPLETED {

  open(C, "Completed_haplotypes") || print "Can't read completed haplotypes list\n";
  while($C = <C>) {
    chomp($C);
    $Completed{$C} = 1
  }
  close C
}



