Here we provide instructions to run the CHI_script.pl script. We also define all variables that need to be set by the user.

Line 8, $window_scan: set to 1 if you choose to analyze msms simulations in sliding windows.
Line 11, @windows: this is an array holding the window boundaries if you choose to analyze simulations via a sliding window method. 
Line 15 $MinCount: This variable describes the minimum number of individuals that a SNP must be present in to be included in the analysis. 
Line 16, $neutral: set to 1 if doing neutral simulations and 0 if you are modeling selection.
Line 17, $RequireComplete: This needs to be set to 1 when only considering complete sweeps, 0 otherwise.
Line 18, $IncSweep: If you are studying partial sweeps, set this variable to 1. The next 2 variables (line 19 and 20) will need to be set as well.
Line 19, $MinAlleleCount: If simulating partial sweeps, the minimum beneficial allele count can be set here. For instance, if simulating a population with 50 sampled individuals and your desired allele frequency range is 46%-54%, this variable would be set to 23.
Line 20, $MaxCount: If simulating partial sweeps, the maximum beneficial allele count can be set here. If simulating a population with 50 sampled individuals and your desired allele frequency range is 46%-54%, this variable would be set to 27.
Line 21, $MinTract: The minimum haplotype sharing required to be included in the total population sharing. In our analyses, we used 10% of the window size (500 bp for high N_e and 10000 bp for low N_e).
Line 22, $XPEHHreps: The number of replicates that you wish to calculate XPEHH on. Simply set to 0 if you do not want to do an XPEHH analysis.
Line 23, $InitialSampleSize1: Here, give the total number of sampled chromosomes from population 1 that you are simulating. See paper for details on why extra chromosomes may need to be simulated for complete sweeps.
Line 24, $FinalSampleSize: The sample size of population 1. When simulating partial sweeps, this variable will be the same as $InitialSampleSize1.
Line 25, $InitialSampleSize2: The number of sampled chromosomes in population 2.
Line 30, $CommandFile: You will need to provide an msms command line. This can be stored in a .txt file, the name of which will be stored in this variable.
Line 31, $InputFile: msms will output the simulation into a file name of your choice. The name will be stored in this variable.
Line 52, $ZipFile: Name zipped file containing msms and xpehh executable programs. The programs will be found in the bin folder of the unzipped file.

To run this script, three files will need to be in the same directory:
1. CHI_script.pl
2. msms command file (1 line .txt file containing the command line for the msms simulation)
3. zipped moms and xpehh file, can be downloaded from the github.

From command line, simply set your working directory to the directory with the three above files and type “perl CHI_script.pl”
