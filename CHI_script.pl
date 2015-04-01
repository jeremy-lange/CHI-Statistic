#!/usr/bin/perl -w
use strict;
use Cwd;

#script to:
# launch msms simulations from a file containing the command line.
# analyze simulated data by calculating CHI, FST, and XPEHH for each replicate (the latter by launching a separate program)
my $window_scan=1;#set to 1 if the calculation of CHI along sliding windows is desired
my @windows; #initialize windows array
if($window_scan==1){
  @windows=([0,0.1],[0.05,0.15],[0.1,0.2],[0.15,0.25],[0.2,0.3],[0.25,0.35],[0.3,0.4],[0.35,0.45],[0.4,0.5],
    [0.45,0.55],[0.5,0.6],[0.55,0.65],[0.6,0.7],[0.65,0.75],[0.7,0.8],[0.75,0.85],[0.8,0.9],[0.85,0.95],[0.9,1]);#array of window start and stop positions
}

my $MinCount = 2;  #if fewer than this number of individuals carry the less common allele at a site, that polymorphism is ignored
my $neutral=0; #set to 1 for a neutral locus, 0 for selection
my $RequireComplete=1; #set to 1 to ensure complete sweeps only, 0 otherwise
my $IncSweep=0; #Set to 1 if you only want to consider incomplete sweeps with final adaptive allele frequencies within a certain range (otherwise set to 0)
my $MinAlleleCount=13;#The minimum number of sampled chromosomes in population 1 that carry the adaptive allele (if it's lower, skip this replicate).  Only matters if $IncSweep = 1.
my $MaxCount=17;#The maximum number of sampled chromosomes in population 1 that carry the adaptive allele (if it's higher, skip this replicate).  Only matters if $IncSweep = 1.
my $MinTract=500;#threshold bp length to analyze identical pairs of haplotypes. In our analyses, 500 for high Ne cases and 10,000 for low Ne cases
my $XPEHHreps = 0;  #Only evaluate XPEHH for this many replicates
my $InitialSampleSize1=52;#For neutral case, just simulate with final sample size (no extras). See paper for reason to simulate extra chromosomes.
my $FinalSampleSize = 50;  #the number of alleles to use in the analysis, for both populations
my $InitialSampleSize2 = 50;  #simulate extra pop2 alleles so we can select those without the adaptive allele
my $window_count;#this as well as next two variables need to be initialized if running a sliding window analysis
my $window_start;
my $window_stop;

my $CommandFile="msms_cmd.txt";#name of one line .txt file containing the msms simulation command line
my $InputFile="msms_output.txt";#name of msms simulation output file

#Declare some stuff
my $i = 0;
my $j = 0;
my $s = 0;
my $TotalPops = 0;
my $TotalSites = 0;
my $reps = 0;
my $median = 0;
my $denominator = 0;
my $NumberToTrim1 = $InitialSampleSize1 - $FinalSampleSize;
my $NumberToTrim2 = $InitialSampleSize2 - $FinalSampleSize;
my @InputAoA = ();
my @InputInd = ();
my @PopSizes = ();
my @LengthSumsA = ();
my @LengthSumsB = ();
my @denominators = ();
my @NonAdaptive = ();

my $ZipFile = "msms_xpehh.zip";#our analyses were run on a computing cluster. We needed to unpack the xpehh and msms programs for each analysis
#unpack the zip file containing msms and a few other files
my $pwd = cwd();#variable storing path of current working directory
system ("unzip $ZipFile");
chdir ("zip");
chdir ("bin");
my $msms_wd = cwd();
system ("chmod +x msms");
system ("chmod +x xpehh");
chdir("$pwd");

#Get msms command line from file, then run msms simulations
open F, "<$CommandFile"
  or die "Can't open input file \"$CommandFile\" ($!)";
my $msmsCommand = (<F>);
close F;
chomp ($msmsCommand);
chdir($msms_wd);
system ($msmsCommand);#this runs the msms simulation

if($window_scan==1){
  for($window_count=0;$window_count<@windows;$window_count++){

$window_start=$windows[$window_count][0];#position along chromosome of the start of the analysis window
$window_stop=$windows[$window_count][1];#position along chromosome of the stop of the analysis window

$reps = 0;

#reset some variables to splice simulations by window
my @InputAoA = ();#multidimensional array containing the simulation (total sample size x # of segregating sites)
my @InputInd = ();
my @PopSizes = ();
my @LengthSumsA = ();
my @LengthSumsB = ();
my @denominators = ();
my @NonAdaptive = ();
my $LengthSum=0;
my @PopComps = ();
##############################################

chdir ("$msms_wd");

#READ THE INPUT DATA

open I, "<$InputFile"
  or die "Can't open input file \"$InputFile\" ($!)";

#get info from command line
$_ = (<I>);
chomp;
my @line = split;
if ($line[0] =~ m/msms/){
  shift @line;
}
my $TotalInds = $line[1];#number of individuals (rows of @InputAoA) that are simulated
my $TotalReps = $line[2];#number of replicates simulated
for ($i = 2; $i < @line; $i++){
  if ($line[$i] =~ m/\-r/){
    $TotalSites = $line[$i+2];#total number of sites simulated
    $TotalSites=($window_stop-$window_start)*$TotalSites;#rescaled by window size
    next;
  }
  if ($line[$i] =~ m/\-I/){
    $TotalPops = $line[$i+1];#number of populations simulated, this program only works with 2
    for ($j = 0; $j < $TotalPops; $j++){
      $s = $i + $j + 2;
      push @PopSizes, $line[$s];#fills @PopSizes with number of individuals simulated in each population
    }
    last;
  }
}
print "From input file $InputFile, expecting data for $TotalInds individuals from $TotalPops populations, for $TotalReps replicates of a locus with $TotalSites sites.\n\nReplicates finished: ";

if (($PopSizes[0] != $InitialSampleSize1) || ($PopSizes[1] != $InitialSampleSize2)){
  die "Found simulated sample sizes $PopSizes[0] and $PopSizes[1].  Expected $InitialSampleSize1 and $FinalSampleSize.\n";
}
$PopSizes[0] = $FinalSampleSize;
$PopSizes[1] = $FinalSampleSize;

if ($TotalPops != 2){
  die "This program only works with two populations.  Found $TotalPops\n";
}
my @PopStarts = ();
my @PopStops = ();
my $start = 0;
my $stop = -1;
push @PopStarts, 0;#the first individual of population 1 is row 0 of @InputAoA
for ($i = 0; $i < @PopSizes; $i++){#loops through array of length 2, each element corresponds to number of individuals in population
  $stop = $stop + $PopSizes[$i];
  push @PopStops, $stop;#array of length 2 holding the index number of the row of the last individual in each popoulation
  last if ($i == (@PopSizes - 1));
  $start = $start + $PopSizes[$i];
  push @PopStarts, $start;#array of length 2 holding the index number of the row of the first individual in each popoulation
}
my $comps = 0;
for ($i = 0; $i < @PopSizes; $i++){
  $comps = 0;
  for ($j = $PopSizes[$i] - 1; $j > 0; $j--){
    $comps += $j;
  }
  push @PopComps, $comps;#array of length 2 containing number of pairwise comparisons in each population
}
my $BtwnComps = $PopSizes[0] * $PopSizes[1];


#Declare some more stuff
my $TargetSite = 0;
my $TargetLoc = 0.50000;#location of the selected site, corresponds to msms "-Sp" flag
my $d = 0;
my $p = 0;
my $count = 0;
my $length = 0;
#my $LengthSum = 0;
my $ratio = 0;
my $pi = 0;
my $diffs = 0;
my $pos = 0;
my $Dxy = 0;
my $FST = 0;
my $XPEHHmax = 0;
my $skip = 1;
my @positions = ();
my @DiffSites = ();
my @CHIASums = ();
my @CHIARatios = ();
my @MedMinRatios = ();
my @TargetFreqs = ();
my @pis = ();
my @FSTs = ();
my @XPEHHs = ();
my @AdaptiveAlleles = ();
my @NumberAdaptiveAlleles = ();

my $first_element;
my $last_element;

#Input loop
scalar (<I>);
while (<I>){
  chomp;
  $_ =~ s/\s+$//;
  next if (m%/%);
  next if (m/segsites/);
  if (m/positions/){
    @positions = split;
    shift @positions;
    if (@positions == 0){
      push @CHIARatios, 1;
    }
    next;
  }
  if ((m/0/) || (m/1/)){
    @line = split //, $_;
    push @InputAoA, [ @line ];

    #Each time a replicate's data matrix is complete (all individuals present), calculate the CHIA ratios...
    if(@InputAoA == $TotalInds){
      
      #Find the advantageous allele and report its frequency in population 
      if ($neutral == 0){
        @AdaptiveAlleles = ();
        $TargetSite = -1;
        for ($s = 0; $s < @positions; $s++){
          if ($positions[$s] == $TargetLoc){
            $TargetSite = $s;
            last;
          }
        }
        if ($TargetSite >= 0){
          $count = 0;
          @NonAdaptive = ();
          for ($i = 0; $i < $InitialSampleSize1; $i++){
            if ($InputAoA[$i][$s] ne '0'){
              $count++;
              push @AdaptiveAlleles, $InputAoA[$i][$s];
              if ($InputAoA[$i][$s] ne '1'){
                $InputAoA[$i][$s] = 1
              }
            }
            else{
              push @NonAdaptive, $i;
            }
          }
          push @TargetFreqs, $count;
        }
        else{
          push @TargetFreqs, 0;
          print "failed to find adaptive allele for a replicate\n";
        }
        #Find the number of unique adaptive alleles in the population
        @AdaptiveAlleles = sort(@AdaptiveAlleles);
        for ($i = 1; $i < @AdaptiveAlleles; $i++){
          if ($AdaptiveAlleles[$i] eq $AdaptiveAlleles[$i-1]){
            splice @AdaptiveAlleles, $i, 1;
            $i--;
          }
        }
        $i=@AdaptiveAlleles;
        push @NumberAdaptiveAlleles, $i;

        #If simulating incomplete sweeps with a specific range of final frequencies, skip replicates that fail this criterion
        if ($IncSweep == 1){
          if (($count > $MaxCount) || ($count < $MinAlleleCount)){
            print "Skipped a replicate because $count individuals had an adaptive allele\n";
            @InputAoA = ();
            pop @TargetFreqs;
            pop @NumberAdaptiveAlleles;
            next;
          }
        }

        #Filter replicates and individuals to obtain completely swept cases in pop 1 (if we're simulating complete sweeps)
        elsif ($count < ($FinalSampleSize * $RequireComplete)){
          @InputAoA = ();
          print "Skipped a replicate because not enough population 1 individuals had adaptive allele: $count\n";
          pop @TargetFreqs;
          pop @NumberAdaptiveAlleles;
          next;
        }
        for ($i = 0; $i < $NumberToTrim1; $i++){
          if ($i <= (@NonAdaptive - 1)){
            $j = $NonAdaptive[$i];
          }
          else{
            $j = $FinalSampleSize;
          }
          splice @InputAoA, $j, 1;
        }
        #Filter replicates and individuals to yield pop 2 without any adaptive alleles
        $skip = 0;
        for ($i = $FinalSampleSize; $i < @InputAoA; $i++){
          if ($InputAoA[$i][$TargetSite] ne '0'){
            splice @InputAoA, $i, 1;
            $i--;
          }
          if ($i == ($FinalSampleSize * 2)){
            for ($j = $i; $j < @InputAoA; $j++){
              splice @InputAoA, $j, 1;
              $j--;
            }
            last;
          }
          elsif ( ($i == (@InputAoA - 1)) && ($i < (($FinalSampleSize * 2) - 1)) ){
            $skip = 1;
            $j = $i - $FinalSampleSize;
            print "Discarded replicate because too few individuals in population 2 lacked adaptive allele: $j\n";
          }
        }
        if ($skip == 1){
          @InputAoA = ();
          pop @TargetFreqs;
          pop @NumberAdaptiveAlleles;
          next;
        }
        if (@InputAoA != ($FinalSampleSize * 2)){
          $i = @InputAoA;
          die "After filtering, found incorrect number of individuals in data matrix ($i)\n";
        }
      }#closes if($neutral==0) statement
      else{
        push @TargetFreqs, 0;
      }

      for($i=0;$i<@positions;$i++){
        if($positions[$i]>=$window_start){
          $first_element=$i;
          last;
        }
      }
      if($window_count>0){
        splice @positions,0,($first_element);
        for($i=0;$i<@InputAoA;$i++){
          splice(@{$InputAoA[$i]},0,($first_element));
        }
      }
      for($i = 0; $i < @positions; $i++){
        if($positions[$i] >= $window_stop){
          $last_element = ($i-1);
          last;
        }else{
          $last_element=@positions;
        }
      }
  
      splice @positions, ($last_element+1);
      for($i=0;$i<@InputAoA;$i++){
        splice(@{$InputAoA[$i]},($last_element+1));
      }
  
      my $pos_track;
      for($pos_track=0;$pos_track<@positions;$pos_track++){
        $positions[$pos_track]=($positions[$pos_track]-$window_start)/($window_stop-$window_start);
      }

      #Evaluate FST
      #calculate pi for each population
      @pis = ();
      for ($p = 0; $p < @PopSizes; $p++){
        $diffs = 0;
        for ($i = $PopStarts[$p]; $i <= ($PopStops[$p] - 1); $i++){
          for ($j = $i + 1; $j <= $PopStops[$p]; $j++){
            for ($s = 0; $s < @{$InputAoA[0]}; $s++){
              if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
                $diffs++;
              }
            }
          }
        }
        $pi = $diffs / ($PopComps[$p] * $TotalSites);
          push @pis, $pi;
      }
      #calculate Dxy and FST for each population pair
      $diffs = 0;
      for ($i = $PopStarts[0]; $i <= $PopStops[0]; $i++){
        for ($j = $PopStarts[1]; $j <= $PopStops[1]; $j++){
          for ($s = 0; $s < @{$InputAoA[0]}; $s++){
            if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
              $diffs++;
            }
          }
        }
      }
      $Dxy = $diffs / ($BtwnComps * $TotalSites);
      if($Dxy>0){
        $FST = 1 - ((($pis[0] + $pis[1]) / 2) / $Dxy);
        push @FSTs, $FST;
        }
      else{
        push @FSTs, -99999;
      }
    
      #Obtain maximum XPEHH by calling an external program
      #create three data files in the necessary format
      if ($reps < $XPEHHreps){
        open A, ">XPEHH_input_pop1.txt";
        for ($i = 0; $i < $PopStops[0]; $i++){
          for ($s = 0; $s < @{$InputAoA[0]}; $s++){
            print A $InputAoA[$i][$s];
            if ($s < (@{$InputAoA[0]} - 1)){
              print A " ";
            }
            else{
              print A "\n";
            }
          }
        }
        close A;
        open B, ">XPEHH_input_pop2.txt";
        for ($i = $PopStarts[1]; $i < $PopStops[1]; $i++){
          for ($s = 0; $s < @{$InputAoA[0]}; $s++){
            print B $InputAoA[$i][$s];
            if ($s < (@{$InputAoA[0]} - 1)){
              print B " ";
            }
            else{
              print B "\n";
            }
          }
        }
        close B;
        open C, ">XPEHH_input_map.txt";
        for ($s = 0; $s < @positions; $s++){
          $pos = $positions[$s] * $TotalSites;
          print C "$pos\t$pos\t$pos\tA\tC\n";
        }
        close C;
        system("./xpehh -m XPEHH_input_map.txt -h XPEHH_input_pop1.txt XPEHH_input_pop2.txt > XPEHH_output.txt");
        open X, "<XPEHH_output.txt" or die "can not open XPEHH output file\n";
        $XPEHHmax = -99999;
        while (<X>){
          chomp;
          @line = split;
          if ($line[4] > $XPEHHmax){
            $XPEHHmax = $line[4];
          }
        }
        close X;
        push @XPEHHs, $XPEHHmax;
      }

      #evalute CHI ratio
      #first remove singletons from InputAoA and positions
      for ($s = 0; $s < @positions; $s++){
        $count = 0;
        for ($i = 0; $i < ($FinalSampleSize * 2); $i++){
          $count += $InputAoA[$i][$s];
        }
        if ( ($count < ($MinCount - 0.5)) || ($count > (@InputAoA - ($MinCount-0.5))) ){
          splice @positions, $s, 1;
          for ($i = 0; $i < @InputAoA; $i++){
            splice @{$InputAoA[$i]}, $s, 1;
          }
          $s--;
        }
      }

      #For each population, cycle through each pair of sequences
      @CHIASums = ();
      for ($p = 0; $p < @PopSizes; $p++){
        $LengthSum = 0;
        for ($i = $PopStarts[$p]; $i < $PopStops[$p]; $i++){
          for ($j = $i + 1; $j <= $PopStops[$p]; $j++){
            #Record sites that differ between this pair of sequences
            @DiffSites = ();
            push @DiffSites, 0;
            for ($s = 0; $s < @positions; $s++){
              if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
                push @DiffSites, $positions[$s];
              }
            }
            push @DiffSites, 1;

            #Identify tracts of identical sequence longer than the threshold
            for ($d = 0; $d < @DiffSites - 1; $d++){
              $length = (($DiffSites[$d+1] - $DiffSites[$d]) * $TotalSites) - 1;
              if ($length >= $MinTract){
                $LengthSum += $length;
              }
            }
          }
        }
        if ($p == 0){
          push @LengthSumsA, $LengthSum;
        }
        else{
          push @LengthSumsB, $LengthSum;
        }
      }#closes ($p = 0; $p < @PopSizes; $p++) loop
      if ($LengthSumsB[-1] > 0){
        $ratio = ($LengthSumsA[-1] / $PopComps[0]) / ($LengthSumsB[-1] / $PopComps[1]);
      }
      else{
        $ratio = ($LengthSumsA[-1] / $PopComps[0]) / (($MinTract / 2) / $PopComps[1]);
      }
      @InputAoA = ();
      $reps++;
      print "replicate $reps, Ben. Allele Freq. $TargetFreqs[-1], CHIA Ratio $ratio.\n";
    }
  }
}

#evaluate simple ratios
for ($i = 0; $i < @LengthSumsA; $i++){
  if ($LengthSumsB[$i] > 0){
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / ($LengthSumsB[$i] / $PopComps[1]);
  }
  else{
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / (($MinTract / 2) / $PopComps[1]);
  }
  push @CHIARatios, $ratio;
}


#evaluate ratios using median denominator as minimum denominator 
for ($i = 0; $i < @LengthSumsB; $i++){
  $denominator = $LengthSumsB[$i] / $PopComps[1];
  push @denominators, $denominator;
}
@denominators = sort { $a <=> $b } @denominators;
$i = @denominators;
if ( ($i % 2) == 1){
  $i = (@denominators - 1) / 2;
  $median = $denominators[$i];
}
else{
  $i = @denominators / 2;
  $median = ($denominators[$i] + $denominators[$i-1]) / 2;
}
if ($median == 0){
  die "median denominator was zero.  try changing your parameters\n";
}
for ($i = 0; $i < @LengthSumsA; $i++){
  $denominator = $LengthSumsB[$i] / $PopComps[1];
  if ($denominator > $median){
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / $denominator;
    push @MedMinRatios, $ratio;
  }
  else{
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / $median;
    push @MedMinRatios, $ratio;
  }
}

#output to file
chdir ("$pwd");
my $file = 'idhap_' . '_start'. $window_start. '_stop'.$window_stop . '_' . $InputFile;
open O, ">$file";
for ($i = 0; $i < @CHIARatios; $i++){
  print O "$TargetFreqs[$i]\t$NumberAdaptiveAlleles[$i]\t$CHIARatios[$i]\t$MedMinRatios[$i]\t$FSTs[$i]\t";
  if ($i < @XPEHHs){
    print O $XPEHHs[$i];
  }
  print O "\n";
}
close O;

}#closes for loop looping through windows
}#closes if($window_scan==1) statement
else{#if not doing a window scan analysis
$reps = 0;
my @InputAoA = ();#multidimensional array containing the simulation (total sample size x # of segregating sites)
my @InputInd = ();
my @PopSizes = ();
my @LengthSumsA = ();
my @LengthSumsB = ();
my @denominators = ();
my @NonAdaptive = ();
my $LengthSum=0;
my @PopComps = ();
##############################################

chdir ("$msms_wd");

#READ THE INPUT DATA

open I, "<$InputFile"
  or die "Can't open input file \"$InputFile\" ($!)";

#get info from command line
$_ = (<I>);
chomp;
my @line = split;
if ($line[0] =~ m/msms/){
  shift @line;
}
my $TotalInds = $line[1];#number of individuals (rows of @InputAoA) that are simulated
my $TotalReps = $line[2];#number of replicates simulated
for ($i = 2; $i < @line; $i++){
  if ($line[$i] =~ m/\-r/){
    $TotalSites = $line[$i+2];#total number of sites simulated
    next;
  }
  if ($line[$i] =~ m/\-I/){
    $TotalPops = $line[$i+1];#number of populations simulated, this program only works with 2
    for ($j = 0; $j < $TotalPops; $j++){
      $s = $i + $j + 2;
      push @PopSizes, $line[$s];#fills @PopSizes with number of individuals simulated in each population
    }
    last;
  }
}
print "From input file $InputFile, expecting data for $TotalInds individuals from $TotalPops populations, for $TotalReps replicates of a locus with $TotalSites sites.\n\nReplicates finished: ";

if (($PopSizes[0] != $InitialSampleSize1) || ($PopSizes[1] != $InitialSampleSize2)){
  die "Found simulated sample sizes $PopSizes[0] and $PopSizes[1].  Expected $InitialSampleSize1 and $FinalSampleSize.\n";
}
$PopSizes[0] = $FinalSampleSize;
$PopSizes[1] = $FinalSampleSize;

if ($TotalPops != 2){
  die "This program only works with two populations.  Found $TotalPops\n";
}
my @PopStarts = ();
my @PopStops = ();
my $start = 0;
my $stop = -1;
push @PopStarts, 0;#the first individual of population 1 is row 0 of @InputAoA
for ($i = 0; $i < @PopSizes; $i++){#loops through array of length 2, each element corresponds to number of individuals in population
  $stop = $stop + $PopSizes[$i];
  push @PopStops, $stop;#array of length 2 holding the index number of the row of the last individual in each popoulation
  last if ($i == (@PopSizes - 1));
  $start = $start + $PopSizes[$i];
  push @PopStarts, $start;#array of length 2 holding the index number of the row of the first individual in each popoulation
}
my $comps = 0;
for ($i = 0; $i < @PopSizes; $i++){
  $comps = 0;
  for ($j = $PopSizes[$i] - 1; $j > 0; $j--){
    $comps += $j;
  }
  push @PopComps, $comps;#array of length 2 containing number of pairwise comparisons in each population
}
my $BtwnComps = $PopSizes[0] * $PopSizes[1];


#Declare some more stuff
my $TargetSite = 0;
my $TargetLoc = 0.50000;#location of the selected site, corresponds to msms "-Sp" flag
my $d = 0;
my $p = 0;
my $count = 0;
my $length = 0;
#my $LengthSum = 0;
my $ratio = 0;
my $pi = 0;
my $diffs = 0;
my $pos = 0;
my $Dxy = 0;
my $FST = 0;
my $XPEHHmax = 0;
my $skip = 1;
my @positions = ();
my @DiffSites = ();
my @CHIASums = ();
my @CHIARatios = ();
my @MedMinRatios = ();
my @TargetFreqs = ();
my @pis = ();
my @FSTs = ();
my @XPEHHs = ();
my @AdaptiveAlleles = ();
my @NumberAdaptiveAlleles = ();

my $first_element;
my $last_element;

#Input loop
scalar (<I>);
while (<I>){
  chomp;
  $_ =~ s/\s+$//;
  next if (m%/%);
  next if (m/segsites/);
  if (m/positions/){
    @positions = split;
    shift @positions;
    if (@positions == 0){
      push @CHIARatios, 1;
    }
    next;
  }
  if ((m/0/) || (m/1/)){
    @line = split //, $_;
    push @InputAoA, [ @line ];

    #Each time a replicate's data matrix is complete (all individuals present), calculate the CHIA ratios...
    if(@InputAoA == $TotalInds){
      
      #Find the advantageous allele and report its frequency in population 
      if ($neutral == 0){
        @AdaptiveAlleles = ();
        $TargetSite = -1;
        for ($s = 0; $s < @positions; $s++){
          if ($positions[$s] == $TargetLoc){
            $TargetSite = $s;
            last;
          }
        }
        if ($TargetSite >= 0){
          $count = 0;
          @NonAdaptive = ();
          for ($i = 0; $i < $InitialSampleSize1; $i++){
            if ($InputAoA[$i][$s] ne '0'){
              $count++;
              push @AdaptiveAlleles, $InputAoA[$i][$s];
              if ($InputAoA[$i][$s] ne '1'){
                $InputAoA[$i][$s] = 1
              }
            }
            else{
              push @NonAdaptive, $i;
            }
          }
          push @TargetFreqs, $count;
        }
        else{
          push @TargetFreqs, 0;
          print "failed to find adaptive allele for a replicate\n";
        }
        #Find the number of unique adaptive alleles in the population
        @AdaptiveAlleles = sort(@AdaptiveAlleles);
        for ($i = 1; $i < @AdaptiveAlleles; $i++){
          if ($AdaptiveAlleles[$i] eq $AdaptiveAlleles[$i-1]){
            splice @AdaptiveAlleles, $i, 1;
            $i--;
          }
        }
        $i=@AdaptiveAlleles;
        push @NumberAdaptiveAlleles, $i;

        #If simulating incomplete sweeps with a specific range of final frequencies, skip replicates that fail this criterion
        if ($IncSweep == 1){
          if (($count > $MaxCount) || ($count < $MinAlleleCount)){
            print "Skipped a replicate because $count individuals had an adaptive allele\n";
            @InputAoA = ();
            pop @TargetFreqs;
            pop @NumberAdaptiveAlleles;
            next;
          }
        }

        #Filter replicates and individuals to obtain completely swept cases in pop 1 (if we're simulating complete sweeps)
        elsif ($count < ($FinalSampleSize * $RequireComplete)){
          @InputAoA = ();
          print "Skipped a replicate because not enough population 1 individuals had adaptive allele: $count\n";
          pop @TargetFreqs;
          pop @NumberAdaptiveAlleles;
          next;
        }
        for ($i = 0; $i < $NumberToTrim1; $i++){
          if ($i <= (@NonAdaptive - 1)){
            $j = $NonAdaptive[$i];
          }
          else{
            $j = $FinalSampleSize;
          }
          splice @InputAoA, $j, 1;
        }
        #Filter replicates and individuals to yield pop 2 without any adaptive alleles
        $skip = 0;
        for ($i = $FinalSampleSize; $i < @InputAoA; $i++){
          if ($InputAoA[$i][$TargetSite] ne '0'){
            splice @InputAoA, $i, 1;
            $i--;
          }
          if ($i == ($FinalSampleSize * 2)){
            for ($j = $i; $j < @InputAoA; $j++){
              splice @InputAoA, $j, 1;
              $j--;
            }
            last;
          }
          elsif ( ($i == (@InputAoA - 1)) && ($i < (($FinalSampleSize * 2) - 1)) ){
            $skip = 1;
            $j = $i - $FinalSampleSize;
            print "Discarded replicate because too few individuals in population 2 lacked adaptive allele: $j\n";
          }
        }
        if ($skip == 1){
          @InputAoA = ();
          pop @TargetFreqs;
          pop @NumberAdaptiveAlleles;
          next;
        }
        if (@InputAoA != ($FinalSampleSize * 2)){
          $i = @InputAoA;
          die "After filtering, found incorrect number of individuals in data matrix ($i)\n";
        }
      }#closes if($neutral==0) statement
      else{
        push @TargetFreqs, 0;
      }
      
      #Evaluate FST
      #calculate pi for each population
      @pis = ();
      for ($p = 0; $p < @PopSizes; $p++){
        $diffs = 0;
        for ($i = $PopStarts[$p]; $i <= ($PopStops[$p] - 1); $i++){
          for ($j = $i + 1; $j <= $PopStops[$p]; $j++){
            for ($s = 0; $s < @{$InputAoA[0]}; $s++){
              if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
                $diffs++;
              }
            }
          }
        }
        $pi = $diffs / ($PopComps[$p] * $TotalSites);
          push @pis, $pi;
      }
      #calculate Dxy and FST for each population pair
      $diffs = 0;
      for ($i = $PopStarts[0]; $i <= $PopStops[0]; $i++){
        for ($j = $PopStarts[1]; $j <= $PopStops[1]; $j++){
          for ($s = 0; $s < @{$InputAoA[0]}; $s++){
            if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
              $diffs++;
            }
          }
        }
      }
      $Dxy = $diffs / ($BtwnComps * $TotalSites);
      if($Dxy>0){
        $FST = 1 - ((($pis[0] + $pis[1]) / 2) / $Dxy);
        push @FSTs, $FST;
        }
      else{
        push @FSTs, -99999;
      }
    
      #Obtain maximum XPEHH by calling an external program
      #create three data files in the necessary format
      if ($reps < $XPEHHreps){
        open A, ">XPEHH_input_pop1.txt";
        for ($i = 0; $i < $PopStops[0]; $i++){
          for ($s = 0; $s < @{$InputAoA[0]}; $s++){
            print A $InputAoA[$i][$s];
            if ($s < (@{$InputAoA[0]} - 1)){
              print A " ";
            }
            else{
              print A "\n";
            }
          }
        }
        close A;
        open B, ">XPEHH_input_pop2.txt";
        for ($i = $PopStarts[1]; $i < $PopStops[1]; $i++){
          for ($s = 0; $s < @{$InputAoA[0]}; $s++){
            print B $InputAoA[$i][$s];
            if ($s < (@{$InputAoA[0]} - 1)){
              print B " ";
            }
            else{
              print B "\n";
            }
          }
        }
        close B;
        open C, ">XPEHH_input_map.txt";
        for ($s = 0; $s < @positions; $s++){
          $pos = $positions[$s] * $TotalSites;
          print C "$pos\t$pos\t$pos\tA\tC\n";
        }
        close C;
        system("./xpehh -m XPEHH_input_map.txt -h XPEHH_input_pop1.txt XPEHH_input_pop2.txt > XPEHH_output.txt");
        open X, "<XPEHH_output.txt" or die "can not open XPEHH output file\n";
        $XPEHHmax = -99999;
        while (<X>){
          chomp;
          @line = split;
          if ($line[4] > $XPEHHmax){
            $XPEHHmax = $line[4];
          }
        }
        close X;
        push @XPEHHs, $XPEHHmax;
      }

      #evalute CHI ratio
      #first remove singletons from InputAoA and positions
      for ($s = 0; $s < @positions; $s++){
        $count = 0;
        for ($i = 0; $i < ($FinalSampleSize * 2); $i++){
          $count += $InputAoA[$i][$s];
        }
        if ( ($count < ($MinCount - 0.5)) || ($count > (@InputAoA - ($MinCount-0.5))) ){
          splice @positions, $s, 1;
          for ($i = 0; $i < @InputAoA; $i++){
            splice @{$InputAoA[$i]}, $s, 1;
          }
          $s--;
        }
      }

      #For each population, cycle through each pair of sequences
      @CHIASums = ();
      for ($p = 0; $p < @PopSizes; $p++){
        $LengthSum = 0;
        for ($i = $PopStarts[$p]; $i < $PopStops[$p]; $i++){
          for ($j = $i + 1; $j <= $PopStops[$p]; $j++){
            #Record sites that differ between this pair of sequences
            @DiffSites = ();
            push @DiffSites, 0;
            for ($s = 0; $s < @positions; $s++){
              if ($InputAoA[$i][$s] != $InputAoA[$j][$s]){
                push @DiffSites, $positions[$s];
              }
            }
            push @DiffSites, 1;

            #Identify tracts of identical sequence longer than the threshold
            for ($d = 0; $d < @DiffSites - 1; $d++){
              $length = (($DiffSites[$d+1] - $DiffSites[$d]) * $TotalSites) - 1;
              if ($length >= $MinTract){
                $LengthSum += $length;
              }
            }
          }
        }
        if ($p == 0){
          push @LengthSumsA, $LengthSum;
        }
        else{
          push @LengthSumsB, $LengthSum;
        }
      }#closes ($p = 0; $p < @PopSizes; $p++) loop
      if ($LengthSumsB[-1] > 0){
        $ratio = ($LengthSumsA[-1] / $PopComps[0]) / ($LengthSumsB[-1] / $PopComps[1]);
      }
      else{
        $ratio = ($LengthSumsA[-1] / $PopComps[0]) / (($MinTract / 2) / $PopComps[1]);
      }
      @InputAoA = ();
      $reps++;
      print "replicate $reps, Ben. Allele Freq. $TargetFreqs[-1], CHIA Ratio $ratio.\n";
    }
  }
}

#evaluate simple ratios
for ($i = 0; $i < @LengthSumsA; $i++){
  if ($LengthSumsB[$i] > 0){
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / ($LengthSumsB[$i] / $PopComps[1]);
  }
  else{
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / (($MinTract / 2) / $PopComps[1]);
  }
  push @CHIARatios, $ratio;
}


#evaluate ratios using median denominator as minimum denominator 
for ($i = 0; $i < @LengthSumsB; $i++){
  $denominator = $LengthSumsB[$i] / $PopComps[1];
  push @denominators, $denominator;
}
@denominators = sort { $a <=> $b } @denominators;
$i = @denominators;
if ( ($i % 2) == 1){
  $i = (@denominators - 1) / 2;
  $median = $denominators[$i];
}
else{
  $i = @denominators / 2;
  $median = ($denominators[$i] + $denominators[$i-1]) / 2;
}
if ($median == 0){
  die "median denominator was zero.  try changing your parameters\n";
}
for ($i = 0; $i < @LengthSumsA; $i++){
  $denominator = $LengthSumsB[$i] / $PopComps[1];
  if ($denominator > $median){
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / $denominator;
    push @MedMinRatios, $ratio;
  }
  else{
    $ratio = ($LengthSumsA[$i] / $PopComps[0]) / $median;
    push @MedMinRatios, $ratio;
  }
}

#output to file
chdir ("$pwd");
my $file = 'idhap_' . '_' . $InputFile;
open O, ">$file";
for ($i = 0; $i < @CHIARatios; $i++){
  print O "$TargetFreqs[$i]\t$NumberAdaptiveAlleles[$i]\t$CHIARatios[$i]\t$MedMinRatios[$i]\t$FSTs[$i]\t";
  if ($i < @XPEHHs){
    print O $XPEHHs[$i];
  }
  print O "\n";
}
close O;

}#closes else statement of if($window_scan==1)

#Clean up by removing bin and lib subdirectories and zip file
system ("rm -rf zip");
chdir("$pwd");
system ("rm $ZipFile");
system ("rm $CommandFile");