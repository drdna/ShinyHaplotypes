package IsItRepeated;

sub REF {

 # Read BLAST report and identify repeated DNA segments

  my $blastFile = $_[0];

  # Keys are sequence IDs. Values are strings, where the character at position
  # i is "0", "1", or "2" if that position has zero hits, one hit, or multiple
  # hits, respectively.

  my $Q;

  my $S;

  my $BLAST;

  my %alignString;
  my $alignStringRef;
  my %QueryCounts = ();
  my %SubjectCounts = ();

  print "Screening BLAST report for repeats: $blastFile\n";

  ($Q, $S, $BLAST) = split(/\./, $blastFile);			        # capture genome identifiers
    
  open(BLAST, "$blastFile") || die "Can't open BLAST file\n";

  print "Query: $Q\tSubject: $S\n";


  # the following examines each BLAST alignment and counts how often each base position in both query and subject occurs in an aligned segment

  while(my $L = <BLAST>) {

    chomp($L);
    next if $L eq '';
    my @BLAST = split(/\t/, $L);

    my($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm) = @BLAST;
    $qid =~ s/.+?(\d+)$/$1/;
    print "$L\n" if $qid eq ''; 
    # Reverse subject start and end if - strand.
    $spos < $send or ($spos, $send) = ($send, $spos);

    # Initialize CountString Hashes
    exists $QueryCounts{$qid} or $QueryCounts{$qid} = "";
    exists $SubjectCounts{$sid} or $SubjectCounts{$sid} = "";

    # Extend QueryCount and SubjectCount strings with zeros if necessary.
    # The 1+ here is because characters in a string are indexed from 0.

    if (length $QueryCounts{$qid} < $qend) {
      $QueryCounts{$qid} .= "0" x ($qend - length $QueryCounts{$qid});
    }

    if (length $SubjectCounts{$sid} < $send) {
      $SubjectCounts{$sid} .= "0" x ($send - length $SubjectCounts{$sid});
    }

    # Replace 0 with 1 (first hit) and 1 with 2 (repeat hit).  The 1+ here
    # is because the range includes both endpoints.
    substr($QueryCounts{$qid}, $qpos, 1 + $qend - $qpos) =~ tr/01/12/;
    substr($SubjectCounts{$sid}, $spos, 1 + $send - $spos) =~ tr/01/12/;
  }


  close BLAST;

  %alignString = %QueryCounts;
  $QueryCounts = undef;
  $SubjectCounts = undef;

  # Re-open BLAST report and determine number of alignments at position; 0 = 0; 1 = 1; 2 = more than 1

  open(BLAST, "$blastFile") || die "Can't open BLAST file\n";

  print "Creating alignment string\n";

  while(my $L = <BLAST>) {

    chomp($L);

    my($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm) = split(/\t/, $L);
    $qid =~ s/.+?(\d+)$/$1/;
    
    # Skip incomplete or blank lines in the BLAST output
    next if $L eq '';

    ## run subroutine that parses the back trace operations portion of BLAST report

    $alignStringRef = ANALYZE_ALIGNMENTS($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm, \%alignString);
  
  }

  close BLAST;
  return($alignStringRef)

}

sub ANALYZE_ALIGNMENTS {

  my($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm, $alignStringRef) = @_;
      
  # nucleotide positions at the start of each alignment need to be adjusted for the first cycle of the analysis

  $qpos --;
  $spos -- if $spos < $send;
  $spos ++ if $spos < $send;

  my @Align = split(/(\d+)/, $AlignSumm);                       # convert BTOP output to list

  shift @Align;

  foreach my $Alignment (@Align) {                              # loop through BTOP list;

    ($qpos, $spos, $alignStringRef) = ANALYZE_ALIGNED_MATCHES($qid, $sid, $qpos, $spos, $send, $Alignment, $alignStringRef) if $Alignment =~ /\d+/;
    ($qpos, $spos, $alignStringRef) = ANALYZE_ALIGNED_MISMATCHES($qid, $sid, $qpos, $spos, $send, $Alignment, $alignStringRef) if $Alignment =~ /[A-Za-z]+/;
  }

  return($alignStringRef)

}


sub ANALYZE_ALIGNED_MATCHES {

# Loop through aligned bases. Increment value at corresponding Align String position by 1 (to a maximum of 2)

  my($qid, $sid, $qpos, $spos, $send, $Alignment, $alignStringRef) = @_;

  %alignString = %$alignStringRef;

  for(my $m = 0; $m <= $Alignment - 1; $m ++) {                                                 # loop through base positions in aligned segment

    # Use 0 if the substring is empty (if this position is past the last match)

    substr($alignString{$qid}, $m -1, 1) =~ tr/1/2/ if substr($SubjectCounts{$sid}, $m - 1, 1) == 2;

    $qpos--;                                                                               # keep track of query base position
    $spos++ if $spos < $send;
    $spos-- if $spos > $send

  }
  return($qpos, $spos, \%alignString)                                       # return results to ANALYZE_ALIGNMENT subroutine

}


sub ANALYZE_ALIGNED_MISMATCHES {

# Loop through mismatched/misaligned bases. Increment value at corresponding Align String position by 1 (to a maximum of 2)

  my ($qid, $sid, $qpos, $spos, $send, $Alignment, $alignStringRef) = @_;

  %alignString = %$alignStringRef;

  my @MutsList = split(//, $Alignment);				# Split mismatch string into a list of characters
  my $NumMutations = @MutsList;					# Count number of characters in list

  for(my $j = 0; $j <= $NumMutations-2; $j += 2) {			# Loop through mismatches, one pair of bases at a time

    my ($qnucl, $snucl) = ($MutsList[$j], $MutsList[$j+1]);		# Read query base and subject base in mismatch

    if(($qnucl =~ /G|A|T|C/) && ($snucl =~ /G|A|T|C/)) {  		# check that both are simple nucleotide substitutions

      # increment base position counters
      $qpos++;
      $spos++ if $spos < $send;
      $spos-- if $spos > $send;
      substr($alignString{$qid}, $qpos -1, 1) =~ tr/1/2/ if substr($SubjectCounts{$sid}, $spos-1, 1) == 2;

    }

    elsif(($qnucl eq "-") && ($snucl =~ /A|G|T|C/))  {		# if query has deletion, no need to report repeat

      # increment base position counters	 					
      $spos++ if $spos < $send;
      $spos-- if $spos > $send;

    }

    elsif(($qnucl =~ /G|A|T|C/) && ($snucl eq "-")) {		# if subject has deletion, report subject repeat on query align string
			
      $qpos++;
      substr($alignString{$qid}, $qpos -1, 1) =~ tr/1/2/ if substr($SubjectCounts{$sid}, $spos - 1, 1) == 2;

    }

  }

  return ($qpos, $spos, \%alignString)				# return results to ANALYZE_ALIGNMENT subroutine

}

1;
