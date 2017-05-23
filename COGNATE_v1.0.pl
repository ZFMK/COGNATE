#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use 5.010;

use FindBin;
use lib "$FindBin::RealBin/"; # find GAL
use File::Path qw(remove_tree);
use Cwd;

use GAL::Annotation;
use GAL::List;

use Bio::DB::Fasta;
use Number::Format;
use List::Util qw(min max sum0);
use Statistics::Basic qw(unbias median variance stddev vector mean handle_missing_values ipres=6);

# newline & tab variables
my $n = "\n";
my $t = "\t";

# number formatting
  my $en = new Number::Format(-thousands_sep   => ',',
                              -decimal_point   => '.');
                                            

#											--- IDEA: Start, stop codons (canonical, count; rest, codon, count)
#											--- IDEA: gene number histogram binned by position on scaffold


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
=head1 COGNATE - COmparative GeNe AnnoTation charactErizator
   

=head1 INTRODUCTION

 This tool takes a genome's .gff (annotation of protein-coding genes) and corresponding .fa file 
 (usually a genome assembly as nucleotide sequences) to produce a 
 concise overview (summary), lists of measured variables, and a component size overview.
 The output includes data and summary statistics for counts, lengths, 
 GC contents, CpGo/e ratio, etc.

 Uses GAL by Barry Moore; written by Jeanne Wilbrandt
 2016-06
 GAL is available here: https://github.com/The-Sequence-Ontology/GAL
    
=head1 SYNOPSIS
 
 perl COGNATE_x.pl --gff /PATH/feature.gff3 --fasta /PATH/genome.fasta [--name OUT_NAME --batch BATCH_NAME] [OPTIONS]
 or
 perl COGNATE_x.pl --input my_COGNATE.input [--batch BATCH_NAME] [OPTIONS]
 
 Paths have to be absolute.
 
 By default, the validity of the given fasta (DNA, only IUPAC characters) and gff
 file (gff3, all in order) is tested, as well as that all seqIDs of the fasta appear
 in field 1 of the gff (seq name). This can be switched off with --no_filecheck (-f).
 
 It is sensible to give OUT_NAME as some sort of species/genome code. If no --name is
 given, a default will be generated (genome_#). Not recommended.
 
 The option --batch allows to write output of multiple analyzed annotations into one file.
 The set of analyzed annotations/assemblies needs a name (BATCH_NAME).
 
=head1 INSTALLATION

 COGNATE requires no installation. GAL libraries have to be in the same directory as the script. 

=head1 OPTIONS

=head2 Required input options

=over 12

=item --input /PATH/COGNATE.input
 
 gives COGNATE a COGNATE.input file, which contains names, gff and fasta files for each to be analyzed genome.
 Mutually exclusive option with --gff and --fasta.
 
=item --fasta /PATH/genome.fasta

 gives COGNATE one fasta file as input. Requires that also --gff is given.
 Mutually exclusive option with --input.
 
=item --gff /PATH/annotation.gff3
 
 gives COGNATE one gff3 file as input. Requires that also --fasta is given.
 Mutually exclusive option with --input.

=back

=head2 Other options

=over 12

=item -n|--name STRING

 name of the genome/species, used in dir/file naming. If not specified, a default will be set (genome_#)
 
=item -b|--batch STRING

 write summary lines to BATCH_NAME output files
 
=item -w|--workingdir
 
 determine working directory

=item -o|--overwrite
 
 overwrite existing COGNATE directories (for OUT_NAME), no questions asked.

=item -f|--no_filecheck

 switch off validity check of input files

=item -l|--transcript_length INT

 evaluate only the longest (max, default), shortest (min) or median (med) transcript. 
 If there are several transcripts fulfilling the criterion, the first matching one is kept.

=item -t|--count_table
 
 switch on the overview table (printed in STDOUT)

=item --print LIST

 tells COGNATE which files to print (default: all).
 LIST has to be a list of numbers, separated by comma, without white-spaces, where
 the numbers are indices of output files (0 - 20).
 Mutually exclusive with --dont_print.


=item --dont_print LIST
 
 tells COGNATE which files not to print.
 LIST has to be a list of numbers, separated by comma, without white-spaces, where
 the numbers are indices of output files (0 - 20).
 Mutually exclusive with --print.

=item -h|--help	
 
 print help

=back

=head2 INPUT

 An input file for COGNATE should follow this format:
 >>
 # # Indicates a comment line, line will be ignored.
 # example_COGNATE.input (file name)
 # NAME		GFF3			FASTA
 species_1	/PATH/annotation_1.gff3	/PATH/genome_1.fa
 genome_b	/PATH/annotation_b.gff3	/home/user/Desktop/genome_b.fa
 <<
 Lines in this file contain three fields, separated by tabs.
 Names should not include white-spaces.
 

=head2 OUTPUT
 
 Of most parameters both mean and median are given. The mean is often inappropriate due to non-normal distribution 
 of data. Thus, for comparisons, the median is recommended.
 
 Find the indices needed for --print / --dont_print in the lists below.
 
 COGNATE will produce a directory (COGNATE_NAME) in the working directory, containing the following files:
 00	COGNATE_NAME_00_analyzed_transcripts.fa (contains protein sequences of all analyzed transcripts)
 01	COGNATE_NAME_01_summary.tsv
 02	COGNATE_NAME_02_scaffold_general.tsv
 03	COGNATE_NAME_03_scaffold_transcripts.tsv
 04	COGNATE_NAME_04_scaffold_CDSs.tsv
 05	COGNATE_NAME_05_scaffold_exons.tsv
 06	COGNATE_NAME_06_scaffold_introns.tsv
 07	COGNATE_NAME_07_transcript_general.tsv
 08	COGNATE_NAME_08_transcript_CDSs.tsv
 09	COGNATE_NAME_09_transcript_exons.tsv
 10	COGNATE_NAME_10_transcript_introns.tsv
 11	COGNATE_NAME_11_CDSs.tsv
 12	COGNATE_NAME_12_exons.tsv
 13	COGNATE_NAME_13_introns.tsv


 ONLY if BATCH_NAME is set (--batch):
 In the working directory, the following files will be produced (not affected by --overwrite):
 14	COGNATE_NAME_14_batch_general.tsv
 15	COGNATE_NAME_15_batch_scaffold-means.tsv
 16	COGNATE_NAME_16_batch_scaffold-medians.tsv
 17	COGNATE_NAME_17_batch_transcript-means.tsv
 18	COGNATE_NAME_18_batch_transcript-medians.tsv
 19	COGNATE_NAME_19_batch_component_sizes.tsv
 20	COGNATE_NAME_20_batch_bash-commands.txt



=head1 VERSIONS 

=head2 1.0
 
 
	 
=cut


my $version = "COGNATE_v1.0.pl";

my $synopsis = "
Usage:
perl $version --gff /PATH/feature.gff3 --fasta /PATH/genome.fasta [--name OUT_NAME --batch BATCH_NAME] [OPTIONS]
or
perl $version --input my_COGNATE.input [--batch BATCH_NAME] [OPTIONS]

If you need more help type:
perl $version --help
";

my $usage = "
############## USAGE ####################

____SYNOPSIS____

perl $version --gff /PATH/feature.gff3 --fasta /PATH/genome.fasta [--name OUT_NAME --batch BATCH_NAME] [OPTIONS]
 or
perl $version --input my_COGNATE.input [--batch BATCH_NAME] [OPTIONS]

 ATTENTION: Make sure that fasta-headers match field1 (location/scaffold) in the gff!
 Paths have to be absolute.

____OPTIONS____
 Mandatory:
   --input COGNATE.input	Input config file, format see below
     or
   --fasta			Fasta file including absolute path for analysis
   --gff			GFF3 file including absolute path for analysis

 Optional:
   -n|--name STRING		Name of the genome/species, used in dir/file naming (OUT_NAME) (ignored when --input is used)
   -w|--workingdir			Determine working directory (Default: current)
   -o|--overwrite			Overwrite existing COGNATE directories (for OUT_NAME), no questions asked.
   -f|--no_filecheck		Switch off validity check of input files
   -b|--batch STRING		Write summary line to BATCH_NAME output files
   --print LIST			Comma-separated list (e.g., 1,2,3) of output file indices (see below), which will be printed
     or
   --dont_print LIST		Comma-separated list (e.g., 1,2,3) of output file indices (see below), which will not be printed
   -l|--transcript_length STRING	Evaluate only the max=longest (default), min=shortest or mid=middle transcript.
   -t|--count_table		Switch on the overview table (printed in STDOUT)
   -h|--help 			Print this (more info: perldoc $version)

 Recommendation: NAME should include a form of species code.

____DESCRIPTION____

COGNATE analyzes a given fasta file and respective gene anotation (gff3) for each genome
(line in COGNATE.input) and outputs data as well as summary statistics for counts, lengths, 
GC contents, CpGo/e ratio, etc.

____INPUT____

An input file for COGNATE should follow this format:
>>
# # Indicates a comment line, line will be ignored.
# example_COGNATE.input (file name)
# NAME	GFF3	FASTA
species_1	/PATH/annotation_1.gff3	/PATH/genome_1.fa
genome_b	/PATH/annotation_b.gff3	/PATH/genome_b.fa
<<
Lines in this file contain three fields, separated by tabs.
Names should not include white-spaces.

____OUTPUT____

The following lists contain the indices for print-control (--print / --dont_print).

The script will produce a directory (COGNATE_NAME) containing the following files:
 00	COGNATE_NAME_00_analyzed_transcripts.fa
 01	COGNATE_NAME_01_summary.tsv
 02	COGNATE_NAME_02_scaffold_general.tsv
 03	COGNATE_NAME_03_scaffold_transcripts.tsv
 04	COGNATE_NAME_04_scaffold_CDSs.tsv
 05	COGNATE_NAME_05_scaffold_exons.tsv
 06	COGNATE_NAME_06_scaffold_introns.tsv
 07	COGNATE_NAME_07_transcript_general.tsv
 08	COGNATE_NAME_08_transcript_CDSs.tsv
 09	COGNATE_NAME_09_transcript_exons.tsv
 10	COGNATE_NAME_10_transcript_introns.tsv
 11	COGNATE_NAME_11_CDSs.tsv
 12	COGNATE_NAME_12_exons.tsv
 13	COGNATE_NAME_13_introns.tsv


ONLY if BATCH_NAME is set (--batch):
In the working directory, the following files will be produced.
These files are not affected by --overwrite.
 14	COGNATE_NAME_14_batch_general.tsv
 15	COGNATE_NAME_15_batch_scaffold-means.tsv
 16	COGNATE_NAME_16_batch_scaffold-medians.tsv
 17	COGNATE_NAME_17_batch_transcript-means.tsv
 18	COGNATE_NAME_18_batch_transcript-medians.tsv
 19	COGNATE_NAME_19_batch_component_sizes.tsv
 20	COGNATE_NAME_20_batch_bash-commands.txt


####################################################
Uses GAL by Barry Moore; written by Jeanne Wilbrandt
2016-06

";


print $n;


# Get options
my $out_name = '';
my ($help, $table_switch, $def_overwrite, $batch_name, $no_filecheck,
	$workingdir, @print, @dont_print, $input, $transcript_choice, $fa_in, $gff_in,
	);

my $opt_success = GetOptions('help|h' 				=> \$help,
							'out_name|name|n=s'		=> \$out_name,
							'fasta|fa|fa_in=s'		=> \$fa_in,
							'gff3|gff|gff_in=s'		=> \$gff_in,
							'overwrite|o' 			=> \$def_overwrite,
							'no_filecheck|f'		=> \$no_filecheck,
							'transcript_length|l=s'	=> \$transcript_choice,
							'batch|b=s' 			=> \$batch_name,
							'count_table|t'			=> \$table_switch,
							'workingdir|wd|w=s'		=> \$workingdir,
							'print=s'				=> \@print,
							'dont_print=s'			=> \@dont_print,
							'input|in=s'			=> \$input,
			      );
if (! $opt_success) {
    die "FATAL: Command_line_parse_error.$n$synopsis";
}

# Print help if needed
if ($help){
	print $usage;
	exit;
}




#### Option checks, output defaults and choices ######################################


# Check for required options

	if (!$fa_in and !$gff_in) {
		die "FATAL: Option --input COGNATE.input is required.$n$synopsis" unless $input; 
	}
	if (!$input) {
		die "FATAL: Option --input COGNATE.input or both --fa FASTA and --gff GFF3 are required.$n$synopsis" unless $fa_in and $gff_in;
	}


# Input file and name
	
	## Input and name are mutually exclusive (names are given in input file)
	die "FATAL: --input FILE and --name STRING are mutually exclusive options. Names should be given in COGNATE.input.$n",
		"Alternatively, give --gff3, --fa and maybe --name.$n$synopsis" if ($input and $out_name);
		
	## Input and fa/gff in are mutually exclusive
	die "FATAL: --input FILE and --fasta FILE + --gff FILE are mutually exclusive options. Choose only one way to submit input files.$n",
		"$synopsis" if ($input and $fa_in or $input and $gff_in);


# ARGV leftovers
	die "FATAL: Something wrong with your input line. \@ARGV contains leftovers:$n - ", 
		join ("$n - ", @ARGV), $n, "Maybe you had spaces in your --print list?$n$synopsis" if @ARGV; 


# Printing

	## @print and @dont_print are mutually excluding print options
	die "FATAL: --print LIST and --dont_print LIST are mutually exclusive options. Define only one of these!$n$synopsis" if (@print and @dont_print);
	
	## Set file print switches (all files are printed by default)
	my @print_file = set_print_files(\@print, \@dont_print);


# Check transcript choice
	$transcript_choice ||= 'max'; # set to default if no value is given
	die "FATAL: Transcript choice must be max (longest), min (shortest) or mid (median)!$n$usage" 
		if ($transcript_choice !~ /max|min|mid/);

# TODO progress?

#### START #############################################################


# Get input files
	
	## Get name (key) and pairs of data (.gff3 and .fa) -> %input = species => (.gff3, .fa)
	my %input;
	
	## If COGNATE.input is used
	if ($input) {
		%input = get_input_data($input);
	}
	
	## If FA and GFF are given
	elsif ($fa_in and $gff_in) {
		my ($feature_file, $fasta_file) = assign_files(($gff_in, $fa_in));
		$input{$out_name} = [$feature_file, $fasta_file];
	}
	
	else {
		die "FATAL: Option --input COGNATE.input or both --fa FASTA and --gff3 GFF3 are required.$n$synopsis"
	}


# Settle working dir

	## Check and change to workingdir
	goto_workingdir($workingdir) if $workingdir;
		
	## Get current dir as variable
	my $cwd = getcwd;



foreach $out_name (keys %input) {
	
		my $now = time;
		
		
		my ($feature_file, $fasta_file) = @{$input{$out_name}};
	
		#### TODO only fasta?!
		
		## Check file validity
		if (!$no_filecheck) {
			my $headers_REF = fasta_validity($fasta_file);
			gff_validity($feature_file, $headers_REF);
		}
	
	
	# Prepare output directory 
	
		## Use standard name when no name is assigned
		$out_name = set_outname($out_name);
	
		## Remove output dir if present and desired
		remove_dir($out_name, $def_overwrite);
	
		## Make output dir and change to it
		goto_outputdir($out_name);
		
	
	# Prepare output files (open file handles)
	
		## Outfile for analyzed transcripts (in aa)
		my $out_fa 				= open_FH('COGNATE_'.$out_name, '00-analyzed_transcripts', 	">", 'fa')	if $print_file[ 0];
	
		## Outfile for summaries (counts, means, etc)
		my $out_summary 		= open_FH('COGNATE_'.$out_name, '01-summary', ">", 'tsv')					if $print_file[ 1];
		
		## Scaffold data ID ref (general, transcripts, CDSs, exons, introns)
		my $OUT_s_general		= open_FH('COGNATE_'.$out_name, '02-scaffold_general', 		">", 'tsv')	if $print_file[ 2];
		my $OUT_s_t				= open_FH('COGNATE_'.$out_name, '03-scaffold_transcripts', 	">", 'tsv')	if $print_file[ 3];
		my $OUT_s_cds			= open_FH('COGNATE_'.$out_name, '04-scaffold_CDSs', 			">", 'tsv')	if $print_file[ 4];
		my $OUT_s_e				= open_FH('COGNATE_'.$out_name, '05-scaffold_exons', 			">", 'tsv')	if $print_file[ 5];
		my $OUT_s_i				= open_FH('COGNATE_'.$out_name, '06-scaffold_introns', 		">", 'tsv')	if $print_file[ 6];
		
		## Transcript data ID ref (general, CDSs, exons, introns)
		my $OUT_t_general		= open_FH('COGNATE_'.$out_name, '07-transcript_general', 		">", 'tsv')	if $print_file[ 7];
		my $OUT_t_cds			= open_FH('COGNATE_'.$out_name, '08-transcript_CDSs', 			">", 'tsv')	if $print_file[ 8];
		my $OUT_t_e				= open_FH('COGNATE_'.$out_name, '09-transcript_exons', 		">", 'tsv')	if $print_file[ 9];
		my $OUT_t_i				= open_FH('COGNATE_'.$out_name, '10-transcript_introns', 		">", 'tsv')	if $print_file[10];
		
		## CDS, exon and intron data, ID ref
		my $OUT_cdss			= open_FH('COGNATE_'.$out_name, '11-CDSs', 	">", 'tsv')					if $print_file[11];
		my $OUT_exons			= open_FH('COGNATE_'.$out_name, '12-exons', 	">", 'tsv')					if $print_file[12];
		my $OUT_introns			= open_FH('COGNATE_'.$out_name, '13-introns', 	">", 'tsv')					if $print_file[13];
	
		## Outfiles for summary lines if batch name is given
		my ($OUT_batch_general, $OUT_batch_s_means, $OUT_batch_s_medians, $OUT_batch_t_means, $OUT_batch_t_medians, $out_compsize, $out_bash);
			
		if ($batch_name) {
			## General data
			$OUT_batch_general		= open_FH('../COGNATE_'.$batch_name, '14-batch_general', 		">>", 'tsv')	if $print_file[14];
			
			## Scaffold meand and medians
			$OUT_batch_s_means		= open_FH('../COGNATE_'.$batch_name, '15-batch_scaffold-means',	">>", 'tsv')	if $print_file[15];
			$OUT_batch_s_medians	= open_FH('../COGNATE_'.$batch_name, '16-batch_scaffold-medians',">>", 'tsv')	if $print_file[16];
		
			## Trnascript means and medians
			$OUT_batch_t_means		= open_FH('../COGNATE_'.$batch_name, '17-batch_transcript-means',	">>", 'tsv')	if $print_file[17];
			$OUT_batch_t_medians	= open_FH('../COGNATE_'.$batch_name, '18-batch_transcript-medians',	">>", 'tsv')	if $print_file[18];
		
			## Component sizes
			$out_compsize 			= open_FH('../COGNATE_'.$batch_name, '19_batch_component_sizes', 	">>",	'tsv')	if $print_file[19];
		
			## Bash commands
			$out_bash 				= open_FH('../COGNATE_'.$batch_name, '20_batch_bash-commands', 		">>", 'txt')	if $print_file[20];
		}
	
	# Load annotation
	
		## Get annotation
		my $annotation	= GAL::Annotation->new($feature_file, $fasta_file);
					      
		## Get Fasta and fasta headers
		my $fasta_db	= Bio::DB::Fasta->new($fasta_file);
		my @scaff_IDs	= $fasta_db->get_all_primary_ids;
	
		## Get annotated features and types
		my $features	= $annotation->features;
		my @types		= $features->types->uniq;
	
	
	# Output count as table (STDOUT) if desired
	if ($table_switch){
		for my $column (qw(seqids sources types strands phases)) {
		  print $features->$column->count_table($column, qw(count));
		} 
	}
	
	
############# GETTING DATA #############################################
	print "# Getting data...$n";
	
	
   
############# SCAFFOLDS and ASSEMBLY #############
	
	# Prepare assembly variables
		my $assembly_length;
		my $assembly_Ns;
		my $assembly_GCs;
		my $assembly_GCSs;
		my $assembly_ambigs;
		my $assembly_nr_C;
		my $assembly_nr_G;
		my $assembly_nr_CG;
		#my $stream  = $fasta_db->get_PrimarySeq_stream;
	
	
	# initiate storage hash for scaffold
		my %scaff_features;
			# %scaff_features = scaff_ID => 
			#	0			scaff_ID,
			#	1			scaff length, 
			#	2			scaff length - N,
			#	3			scaff GC,
			#	4			scaff GC noAm,
			#	5			count transcripts, 
			#	6			added transcript length,
			#	7			transcript coverage,
			#	8			transcript_density, 
			#	9			count cdss,
			#	10			added cds length, 
			#	11			cds coverage,
			#	12			cds density,
			#	13			count exons,
			#	14			added exon length, 
			#	15			exon coverage,
			#	16	 		exon density,
			#	17			count introns,
			#	18			added intron length,
			#	19			intron coverage,
			#	20			intron density
			#	21			transcript_plusstrand
			#	22			transcript_minusstrand
			#	23			transcript_strand_ratio
	
	
	# Get scaffold parameters
		foreach my $scaff_ID (@scaff_IDs) {
			
			## Get length from fasta
			my $scaff_length   		= $fasta_db->length($scaff_ID);
			
			## Calculate length without Ns
			my $sequence 			= $fasta_db->seq($scaff_ID);
			my $Ns					= $sequence =~ tr/Nn/Nn/;
			my $scaff_length_noN	= $scaff_length - $Ns;
			
			## Calculate scaffold GC and GC without ambiguity
				### Get number of bases
				my $scaff_ambigs 		= $sequence =~ tr/NRYKMBDHVnrykmbdhv/NRYKMBDHVnrykmbdhv/;
				my $scaff_GCs			= $sequence =~ tr/GCgc/GCgc/;
				my $scaff_GCSs			= $sequence =~ tr/GCSgcs/GCSgcs/;
				
				### Calculate ratios (=contents)
				my $scaff_GCratio		= sprintf("%.2f" , $scaff_GCs/$scaff_length*100);
				my $scaff_GCnoAmratio 	= 'NA'; 
				if ($scaff_length-$scaff_ambigs > 0) {
					my $scaff_GCnoAmratio	= sprintf("%.2f" , $scaff_GCSs/($scaff_length-$scaff_ambigs)*100);
				}
			
			## Count C, G and CG (di)nucleotides for CpG o/e
				my $nr_C 		= $sequence =~ s/C/C/gi;
				my $nr_G 		= $sequence =~ s/G/G/gi;
				my $nr_CG		= $sequence =~ s/CG/CG/gi;
				
			## Store ID and lengths (and initiate array elements for rest of data)
			push @{$scaff_features{$scaff_ID}}, $scaff_ID, $scaff_length, $scaff_length_noN, $scaff_GCratio, $scaff_GCnoAmratio, 
												0, 0, 0, 0, 0, 0, 0, 0,
												0, 0, 0, 0, 0, 0, 0, 0,
												0, 0, 0,
												;
			
			## Sum up values for assembly parameters
			$assembly_length 	+= $scaff_length;
			$assembly_Ns		+= $Ns;
			$assembly_ambigs	+= $scaff_ambigs;
			$assembly_GCs		+= $scaff_GCs;
			$assembly_GCSs		+= $scaff_GCSs;
			$assembly_nr_C		+= $nr_C;
			$assembly_nr_G		+= $nr_G;
			$assembly_nr_CG		+= $nr_CG;
			
		}
	
		### Further scaffold data will be collected in the transcript section... 


	
	# Calculate assembly parameters
		## Total assembly length: $assembly_length
	
		## Calculate length without Ns
		my $assembly_length_noN = $assembly_length - $assembly_Ns;
		
		## Calculate GC ratios (=contents)
		my $assembly_GCratio 		= sprintf("%.2f" , $assembly_GCs/$assembly_length*100);
		my $assembly_GCnoAmratio 	= sprintf("%.2f" , $assembly_GCSs/($assembly_length - $assembly_ambigs)*100);
	
		## Assembly CpG o/e
		
			### Calculate frequencies of C, G and CG
				my $freq_C 		= $assembly_nr_C/$assembly_length;
				my $freq_G 		= $assembly_nr_G/$assembly_length;
				my $freq_CG		= $assembly_nr_CG/($assembly_length - 1);
				
				my $freq_C_noN	= $assembly_nr_C/$assembly_length_noN;
				my $freq_G_noN	= $assembly_nr_G/$assembly_length_noN;
				my $freq_CG_noN	= $assembly_nr_CG/($assembly_length_noN - 1);
				
				
				my $assembly_CpGoe 		= sprintf "%.3f", ($freq_CG / ($freq_C * $freq_G));
				my $assembly_CpGoe_noN 	= sprintf "%.3f", ($freq_CG_noN / ($freq_C_noN * $freq_G_noN));
	
	
################ TRANSCRIPTS #######################
	print "# Iterating over transcripts... (CDSs, exons, introns)$n";
	
=head1 NOTES
	
=head2 Transcript features
	 
 The term 'transcript' instead of 'mRNA' is used in descriptions, as it is more 
 inclusive (exons, CDSs, introns still contained). 
 Technically, mRNAs are annotated if UTR annotation is given.
	 
 This tool only evaluates one transcript per gene, thus no isoforms although 
 they are given in NCBI gffs. The user may choose whether the 
 longest (= max, default), shortest (= min) or median (=mid) transcript 
 shall be analyzed. 
	 
=head2 Protein length

 Protein length is calculated without stop codon (*).

=cut
	
	# Prepare variables
	
		my %transcript_features; 	
			# %transcript_features = id => 
			# 	0			transcript_ID,
			#	1	general		transcript_length, 
			#	2			transcript_gc, 
			#	3			transcript_gc_noam, 
			#	4			transcript_CpGoe,
			#	5			transcript_strand,
			#	6			protein_length, 
			#	7			altsplice_count
			#	8 	CDSs		cds_count, 
			#	9			added_cds_length, 
			#	10			med_cds_length, 
			#	11			avg_cds_length, 
			#	12			med_cds_gc, 
			#	13			avg_cds_gc, 
			#	14			med_cds_gc_noam, 
			#	15			avg_cds_gc_noam, 
			#	16			med_cds_CpG_oe,
			#	17			cds_coverage, 
			#	18			cds_density,
			#	19	 Exons		e_count,
			#	20			added e_length,
			#	21			med_e_length,
			#	22			avg_e_length,
			#	23			med_e_gc,
			#	24			avg_e_gc,
			#	25			med_e_gc_noam,
			#	26			avg_e_gc_noam,
			#	27			med_e_CpG_oe,
			#	28			e_coverage,
			#	29			e_density,
			#	30	 Introns		i_count,
			#	31			added i_length,
			#	32			med_i_length,
			#	33			avg_i_length,
			#	34			med_i_gc,
			#	35			avg_i_gc,
			#	36			med_i_gc_noam,
			#	37			avg_i_gc_noam,
			#	38			med_i_CpG_oe,
			#	39			i_coverage,
			#	40			i_density

	
								#			0			1		2	3		4		5
		my %cds_features;		# id => transcript_id, length, gc, gc_noam, CpGoe, strand
		my %exon_features; 		# id => transcript_ID, length, gc, gc_noam, CpGoe, strand
		my %intron_features;	# id => transcript_ID, length, gc, gc_noam, CpGoe, strand
		
		# Store ID and corresponding scaffold ID of all analyzed transcripts
		my %scaff_of_transcript; # transcript_ID => scaffold ID
		
		# Store locus (scaff_ID:start-stop) for each transcript --> overlap check
		my %transcript_loci;
	
		my $single_cds_genes		= 0;
		my $single_exon_genes		= 0;
		
		my $strand_mix_genes		= 0;
	
		my $non_iso_transcript_nr	= 0;
	
		# Prepare handling of missing scaffolds
		my $scaff_miss_flag 		= 0;
		my %missing_scaffs; # scaff id => @transcript_IDs
			# if scaffold of transcript is missing in fasta: 
			#	- raise miss-flag --> print missing scaffold ids later
			#	- collect scaffold id => @skipped_transcripts
			# 	- go to next transcript (skip)
			# 	- wrap up: print missing scaffold ids
	
	
		# Introns need handmade ids!
		my $intron_id 				= 0; 
	
		# Flags for presence of CDS / exons in complete annotation
		my $overall_cds_presence = $features->search({type => 'CDS'});
		my $overall_exon_presence = $features->search({type=> 'exon'});
	
	
	
	# Get an iterator for all genes
		my $genes = $features->search({type => 'gene'});
		
		## Die if no genes are in the gff, something is wrong/file cannot be used
		die "FATAL ERROR: No genes available in your annotation file! $n" if $genes == 0;
	
	
	
##### Go through all GENES #############
	while (my $gene = $genes->next) {
		
		# Get one transcript
			
			## Get all transcripts assigned to one gene
			my $transcripts = $gene->mRNAs;
			
			my $altsplice_count = $transcripts - 1;
			
			## If there are none, go to next gene
			next if $transcripts == 0;
		
			## Determine which transcript to keep for analysis (only one per gene, no isoforms)
			## Only take one (the longest = max, shortest = min, median = mid) transcript per gene
			my $transcript = choose_transcript($transcript_choice, $transcripts);
		
			## Skip transcript/gene if no seq is available
			next if !$transcript->protein_seq;

		
		# Count number of analyzed transcripts (1 per gene = number of genes)
		++ $non_iso_transcript_nr;
		
		
		# Get transcript parameters
		
			## Get ID (in case none is annotated, make one)
			my $transcript_id;
			my $id_suffix = 'a';
			if ($transcript->feature_id) {
				$transcript_id = $transcript->feature_id;
			}
			else {
				$transcript_id = 'transcript_'.$id_suffix;
				++$id_suffix;
			}
		
			## Get length
			my $transcript_length = $transcript->genomic_length or warn "FAIL: Cannot retrieve genomic length for transcript '$transcript_id'.$n";
			
						
			## Skip transcript/gene if nuc sequence is only N (raise warning)-------------------------------------------TODO as sub?
			my $transcript_seq 		= $transcript->seq;
			my $N_in_seq 			= $transcript_seq =~ tr/Nn/Nn/;
			my $no_N_length			= $transcript_length - $N_in_seq;
			if ($no_N_length == 0) {
				warn "WARN: Transcript '$transcript_id' consists only of Ns!$n";
				next;
			}
			
			## Get strand (prepare strand variables)
			my $transcript_strand = $transcript->strand;
			my $strand_mix_flag = 0;
			
			my $cds_strands; # +/-/! (! if conflicts appeared)
			my $exon_strands; 
			my $intron_strands;
			
			## Get locus for later overlap check
			my $locus = $transcript->locus;
			push @{$transcript_loci{$transcript_id}}, $locus, $transcript_strand;
			
		
=head2 GC content

 For GC content, two types are calculated: total (gc/length) and (non-)ambiguity (gcs/length-NRYKMBDHV) [noAm].
 The latter GC content is not dependent on assembly quality. Both types include softmasked sequences.

=cut
			## Determine GC contents
			my $transcript_GC 		= sprintf("%.2f" , ($transcript->gc_content)*100);
			my $transcript_GCnoAm	= sprintf("%.2f" , ($transcript->gc_content_noAm)*100);
			my $transcript_CpGoe	= calculate_CpG_oe($transcript->seq);
	  
			## Get protein seq
			my $protein 		= $transcript->protein_seq;
			my $protein_length	= length $protein;
			
			## Print protein seq to analyzed_transcripts.fa
			print $out_fa ">${out_name}|${out_name}_$transcript_id\n$protein$n" if $print_file[0];
	
		# Put transcript parameters into feature hash
		push @{$transcript_features{$transcript_id}}, ($transcript_id, $transcript_length, $transcript_GC, 
														$transcript_GCnoAm, $transcript_CpGoe, $transcript_strand, 
														$protein_length, $altsplice_count);
	    
		
		# Scaffold related features/parameters:
			
			## Get ID of scaffold transcript is located on
			my $transcript_scaff = $transcript->seqid;
			$scaff_of_transcript{$transcript_id} = $transcript_scaff;
		
			## Check whether scaffold is in fasta!
			if (grep {$_ eq $transcript_scaff} @scaff_IDs){
				# if scaffold bearing transcript is element of fasta-scaffolds, all is well				
			}
			else { # if not present 
				$scaff_miss_flag = 1;							 # set flag to trigger printing later!
				push @{$missing_scaffs{$transcript_scaff}}, ($transcript_id);  # remember missing id and transcript IDs that were skipped
				next;
			}
		
			## Count transcripts per scaffold
			++$scaff_features{$transcript_scaff}[5];
			
			## Collect overall transcript length per scaffold
			$scaff_features{$transcript_scaff}[6] += $transcript_length;
			
			## Count transcripts on strands per scaffold
			++$scaff_features{$transcript_scaff}[21] if $transcript_strand eq '+';
			++$scaff_features{$transcript_scaff}[22] if $transcript_strand eq '-';
			
		
	
		# Get counts (know whether features are present and count genes with only one CDS/exon)
		my $cds_count		= $transcript->CDSs->all;
		my $cds_presence	= $cds_count;
		++$single_cds_genes if ($cds_count == 1);
		
		my $e_count 		= $transcript->exons->all;
		my $exon_presence 	= $e_count;
		++$single_exon_genes if ($e_count == 1);
		
		my $i_count 	= $transcript->introns->all;
	  
			# add counts to scaffold data
			$scaff_features{$transcript_scaff}[ 9] += $cds_count;
			$scaff_features{$transcript_scaff}[13] += $e_count;
			$scaff_features{$transcript_scaff}[17] += $i_count;
		
	
#### Go through CDSs per transcript #######################
		
		# Prepare variables
		my (@cds_lengths, $added_cds_length, @cds_gcs, @cds_gcs_noam, @cds_CpGoes);
		
		
		# Get iterator for CDSs of the current transcript
		my $cdss = $transcript->CDSs;
	  
		# Get data of CDS if CDSs are present
		if ($cds_presence) {
		
			# Prepare suffixes for hand-made IDs
			my $id_suffix	= 'a';
			my $sec_suffix	= 'a';
		
			# Iterate over each CDS
			while (my $cds = $cdss->next) {
				## Get ID
				my ($cds_id, $id_suffix) = get_or_make_ID($transcript_id, $cds, 'CDS', $id_suffix);
			
				## Get lengths and GC contents
				my ($cds_length, $cds_gc, $cds_gc_noam)	= get_lengths_and_GCs($cds);		
				$added_cds_length +=$cds_length;
				
				my $cds_CpGoe = calculate_CpG_oe($cds->seq);
			
				my $cds_strand = $cds->strand;
				
				## Strandedness (same as transcript? If not: !)
				if ($transcript_strand ne $cds_strand) {
					$cds_strand .= '!';
					$strand_mix_flag = 1;
				}
				
				
				## Store data for later use
				if (exists $cds_features{$cds_id}) {
					$cds_id = $cds_id."_$sec_suffix";
					++$sec_suffix;
				}
			
				push @{$cds_features{$cds_id}}, ($transcript_id, $cds_length, $cds_gc, $cds_gc_noam, $cds_CpGoe, $cds_strand);
				push @cds_lengths, $cds_length;
				push @cds_gcs, $cds_gc;
				push @cds_gcs_noam, $cds_gc_noam;
				push @cds_CpGoes, $cds_CpGoe;
			}
		
			# Collect CDS features
		
				## Calculate medians and averages for CDSs of the transcript
				my $med_cds_length 	= sprintf("%.2f" , my_median(\@cds_lengths)); 
				my $avg_cds_length 	= sprintf("%.2f" , my_mean(\@cds_lengths));
				my $med_cds_gc 		= sprintf("%.2f" , my_median(\@cds_gcs));
				my $avg_cds_gc 		= sprintf("%.2f" , my_mean(\@cds_gcs));
				my $med_cds_gc_noam	= sprintf("%.2f" , my_median(\@cds_gcs_noam));
				my $avg_cds_gc_noam	= sprintf("%.2f" , my_mean(\@cds_gcs_noam));
				my $med_cds_CpGoe	= my_median(\@cds_CpGoes);
				
				## CDS coverage (cds length / transcript length)
				my $cds_coverage =  sprintf("%.15f" , ($added_cds_length / $transcript_length));
				
				## CDS density (cds count / transcript length)
				my $cds_density  = sprintf("%.15f" , ($cds_count / $transcript_length));
				

		  
			# Push into feature hashes
			push @{$transcript_features{$transcript_id}}, ($cds_count, $added_cds_length, $med_cds_length, $avg_cds_length, 
												$med_cds_gc, $avg_cds_gc, $med_cds_gc_noam, $avg_cds_gc_noam, $med_cds_CpGoe, 
												$cds_coverage, $cds_density); 
			$scaff_features{$transcript_scaff}[9] += $added_cds_length;	 
		}
	  
	  
		# If no CDS is present for this transcript, give undef (-> empty cell) to storage hash
		else {
			push @{$transcript_features{$transcript_id}}, (0, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'); 
			$scaff_features{$transcript_scaff}[9] += 0;
		}
	  
	  
	  
###  Go through Exons per transcript ####################
		
		# Get iterator for all exons of the considered transcript
		my $exons = $transcript->exons;
	  
		# If exons are present...
		if ($exon_presence) {
			
			# Prepare variables
			my (@e_lengths, $added_e_length, @e_gcs, @e_gcs_noam, @e_CpGoes);
		
			# Prepare suffixes for hand-made IDs
			my $id_suffix = 'a';
			my $sec_suffix	= 'a';
			
			# Iterate over each exon
			while (my $exon = $exons->next) {
				
				## Get ID
				my ($exon_id, $id_suffix) = get_or_make_ID($transcript_id, $exon, 'Exon', $id_suffix);
				
				## Get lengths and GC contents
				my ($e_length, $e_gc, $e_gc_noam) = get_lengths_and_GCs($exon);
				$added_e_length += $e_length;
				my $e_CpGoe = calculate_CpG_oe($exon->seq);
				
				my $e_strand = $exon->strand;
				## Strandedness (same as transcript? If not: !)
				if ($transcript_strand ne $e_strand) {
					$e_strand .= '!';
					$strand_mix_flag = 1;
				}
			
				## Put into storage for later use
				if (exists $exon_features{$exon_id}) {
					$exon_id = $exon_id."_$sec_suffix";
					++$sec_suffix;
				}
				push @{$exon_features{$exon_id}}, ($transcript_id, $e_length, $e_gc, $e_gc_noam, $e_CpGoe, $e_strand);	
				push @e_lengths, $e_length;
				push @e_gcs, $e_gc;
				push @e_gcs_noam, $e_gc_noam;
				push @e_CpGoes, $e_CpGoe;
			}
			
			
			# Collect exon features
			
				## Calculate medians and averages for exons of the current transcript
				my $med_e_length 	= sprintf("%.2f" , my_median(\@e_lengths));
				my $avg_e_length 	= sprintf("%.2f" , my_mean(\@e_lengths));
				my $med_e_gc 		= sprintf("%.2f" , my_median(\@e_gcs));
				my $avg_e_gc 		= sprintf("%.2f" , my_mean(\@e_gcs));
				my $med_e_gc_noam 	= sprintf("%.2f" , my_median(\@e_gcs_noam));
				my $avg_e_gc_noam 	= sprintf("%.2f" , my_mean(\@e_gcs_noam));
				my $med_e_CpGoe		= my_mean(\@e_CpGoes);
		
				## Exon coverage (exon length / transcript length)
				my $e_coverage =  sprintf("%.15f" , ($added_e_length / $transcript_length));
			
				## Exon density (exon count /transcript length)
				my $e_density  = sprintf("%.15f" , ($e_count / $transcript_length));
				
			
			# Push into feature hash
			push @{$transcript_features{$transcript_id}}, ($e_count, $added_e_length, $med_e_length, $avg_e_length,
											$med_e_gc, $avg_e_gc, $med_e_gc_noam, $avg_e_gc_noam, $med_e_CpGoe,
											$e_coverage, $e_density);
			$scaff_features{$transcript_scaff}[20] += $added_e_length;
		}
		
		# If no exons are present for this transcript, give undef (-> empty cell) to storage hash
		else {
			push @{$transcript_features{$transcript_id}}, (0, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA');
			$scaff_features{$transcript_scaff}[20] += 0;
		}
	  
###  Go through Introns per transcript ####################
	
		# Prepare variables
		my (@i_lengths, $added_i_length, @i_gcs, @i_gcs_noam, @i_CpGoes) = 0;
		
		# Prepare id_suffix for hand-made intron IDs
		my $i_id_suffix = 'a';
		
		
		# Get iterator for introns
		my $introns = $transcript->introns;
	
		# Go through all introns
		while (my $intron = $introns->next) {
			
			# Get lengths and GC contents
			my ($i_length, $i_gc, $i_gc_noam) = get_lengths_and_GCs($intron);
			$added_i_length += $i_length;
			my $i_CpGoe = calculate_CpG_oe($intron->seq);
			
			my $i_strand = $intron->strand;
			## Strandedness (same as transcript? If not: add !)
			if ($transcript_strand ne $i_strand) {
				$i_strand .= '!';
				$strand_mix_flag = 1;
			}
			
			# Put into storage
			if (exists $intron_features{$intron_id}) {
				$intron_id = $intron_id."_$i_id_suffix";
				++$i_id_suffix;
			}
			push @{$intron_features{$intron_id}}, ($transcript_id, $i_length, $i_gc, $i_gc_noam, $i_CpGoe, $i_strand);
			push @i_lengths, $i_length;
			push @i_gcs, $i_gc;
			push @i_gcs_noam, $i_gc_noam;
			push @i_CpGoes, $i_CpGoe;
	    
			# Raise ID for next intron in this transcript!
			++$intron_id;
		}
		
		# Only if introns are present
		if ($i_count != 0) { 
	
			# Collect intron features
			
				## Calculate medians and averages
				my $med_i_length 	= sprintf("%.2f" , my_median(\@i_lengths));
				my $avg_i_length 	= sprintf("%.2f" , my_mean(\@i_lengths));
				my $med_i_gc 		= sprintf("%.2f" , my_median(\@i_gcs));
				my $avg_i_gc 		= sprintf("%.2f" , my_mean(\@i_gcs));
				my $med_i_gc_noam	= sprintf("%.2f" , my_median(\@i_gcs_noam));
				my $avg_i_gc_noam	= sprintf("%.2f" , my_mean(\@i_gcs_noam));
				my $med_i_CpGoe		= my_mean(\@i_CpGoes);
		  
				## Intron coverage (intron length / transcript length)
				my $i_coverage 	= sprintf("%.15f" , ($added_i_length/$transcript_length));
			
				## Intron density (intron count /transcript length)
				my $i_density  = sprintf("%.15f" , ($i_count / $transcript_length));
		  
			# Push into feature hashes
			push @{$transcript_features{$transcript_id}}, ($i_count, $added_i_length, $med_i_length, $avg_i_length, 
												$med_i_gc, $avg_i_gc, $med_i_gc_noam, $avg_i_gc_noam, $med_i_CpGoe,
												$i_coverage, $i_density);
			$scaff_features{$transcript_scaff}[31] += $added_i_length;  
		}
		
		# If no introns are present for this transcript, give undef (-> empty cell) to storage hash
		else { 
			push @{$transcript_features{$transcript_id}}, (0, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA');
			$scaff_features{$transcript_scaff}[31] += 0;
		}
		
		
		
		# Count genes with strand-mix (i.e., where CDS/exons/introns do have differing strands from the transcript)
		if ($strand_mix_flag == 1) {
			++ $strand_mix_genes
		}
	}
	
	
	
########### WRAP UP AND PRINTING #######################################
	print "${t}Found $strand_mix_genes transcripts with mixed strandedness of substructures. Check files 11/12/13 for '!'.$n"
			if $strand_mix_genes > 0;

### SCAFFOLDS ###########################################
			# %scaff_features = scaff_ID => 
			#	0			scaff_ID,
			#	1			scaff length, 
			#	2			scaff length - N,
			#	3			scaff GC,
			#	4			scaff GC noAm,
			#	5			count transcripts, 
			#	6			added transcript length,
			#	7			transcript coverage,
			#	8			transcript_density, 
			#	9			count cdss,
			#	10			added cds length, 
			#	11			cds coverage,
			#	12			cds density,
			#	13			count exons,
			#	14			added exon length, 
			#	15			exon coverage,
			#	16	 		exon density,
			#	17			count introns,
			#	18			added intron length,
			#	19			intron coverage,
			#	20			intron density
			#	21			transcript_plusstrand
			#	22			transcript_minusstrand
			#	23			transcript_strand_ratio
		
	
		# Go through all scaffolds and collect remaining data
			## Coverage	= feat length	/ scaff length
			## Density	= feat number	/ scaff lenght
			## transcript_strand_ratio = +transcript/transcript count *100 : -transcripts/transcript count *100
			
			## Calculated values are directly stored in the storage hash
			
			foreach my $scaff (@scaff_IDs) {
				## Transcript coverage and density
				$scaff_features{$scaff}[7] = sprintf("%.15f" , ($scaff_features{$scaff}[6]/$scaff_features{$scaff}[1]));
				$scaff_features{$scaff}[8] = sprintf("%.15f" , ($scaff_features{$scaff}[5]/$scaff_features{$scaff}[1]));
				
				## CDS coverage and density
				$scaff_features{$scaff}[11] = sprintf("%.15f" , ($scaff_features{$scaff}[10]/$scaff_features{$scaff}[1])) if $overall_cds_presence;
				$scaff_features{$scaff}[12] = sprintf("%.15f" , ($scaff_features{$scaff}[9]/$scaff_features{$scaff}[1])) if $overall_cds_presence;
				
				## Exon coverage and density
				$scaff_features{$scaff}[15] = sprintf("%.15f" , ($scaff_features{$scaff}[14]/$scaff_features{$scaff}[1])) if $overall_exon_presence;
				$scaff_features{$scaff}[16] = sprintf("%.15f" , ($scaff_features{$scaff}[13]/$scaff_features{$scaff}[1])) if $overall_exon_presence;
				
				## Intron coverage and density
				$scaff_features{$scaff}[19] = sprintf("%.15f" , ($scaff_features{$scaff}[18]/$scaff_features{$scaff}[1]));
				$scaff_features{$scaff}[20] = sprintf("%.15f" , ($scaff_features{$scaff}[17]/$scaff_features{$scaff}[1]));
				
				## Ratio of transcript strands
				my $plus_ratio	= 'NA';
				my $minus_ratio = 'NA';
					# if there are transcripts (otherwise division by zero), calculate ratio
				if ($scaff_features{$scaff}[5] > 0) {
					$plus_ratio		= sprintf("%.2f", ($scaff_features{$scaff}[21]/$scaff_features{$scaff}[5])*100);
					$minus_ratio	= sprintf("%.2f", ($scaff_features{$scaff}[22]/$scaff_features{$scaff}[5])*100);
				}
				$scaff_features{$scaff}[23] = $plus_ratio.' : '.$minus_ratio;
			}
	
	
	
###### GET ELEMENTS for PRINTING and final calculations ############
	
	# Elements of scaffolds
		my @scaff_rows_IDsort		= @scaff_features{sort {$a cmp $b} keys %scaff_features};
		my @scaff_columns			= transpose(\@scaff_rows_IDsort);
		#my @scaff_rows_ValSort		= transpose_sort_transpose(\@scaff_rows_IDsort);
	
	# Elements of transcripts
		my @transcript_rows_IDsort 	= @transcript_features{sort {$a cmp $b} keys %transcript_features};
		my @transcript_columns		= transpose(\@transcript_rows_IDsort);
		#
	
	# Elements of CDSs
		my @cds_rows_IDsort			= @cds_features{sort {$a cmp $b} keys %cds_features} 	if $overall_cds_presence;
		my @cds_columns				= transpose(\@cds_rows_IDsort)							if $overall_cds_presence;
		#
	
	# Elements of exons
		my @exon_rows_IDsort		= @exon_features{sort {$a cmp $b} keys %exon_features}	if $overall_exon_presence;
		my @exon_columns			= transpose(\@exon_rows_IDsort)							if $overall_exon_presence;
		#
	
	# Elements of introns
		my @intron_rows_IDsort		= @intron_features{sort {$a cmp $b} keys %intron_features};
		my @intron_columns			= transpose(\@intron_rows_IDsort);
		#
	
	
### ASSEMBLY DATA ############################################
	
		## N50, 75, 90 / L50, 75, 90
			my @sorted_scaff_lengths = sort {$b <=> $a} @{$scaff_columns[1]};
					
			my ($N50, $L50) = assembly_NL(\@sorted_scaff_lengths, $assembly_length*0.5);
			my ($N75, $L75) = assembly_NL(\@sorted_scaff_lengths, $assembly_length*0.75);
			my ($N90, $L90) = assembly_NL(\@sorted_scaff_lengths, $assembly_length*0.9);
	
		## L90(genes)
			
			### Get counts of genes on scaffolds
			my %genes_on_scaff;
			
			foreach my $transcript_ID (keys %scaff_of_transcript) {
				my $scaff_ID = $scaff_of_transcript{$transcript_ID};
				++$genes_on_scaff{$scaff_ID};
			}
			
			### Prepare variables
			my $L90_genes; 		# Counted scaffolds
			my $gene_count;		# Counted genes
			
			### Go through sorted list of scaffold IDs (sorted by length, longest first)
				foreach my $scaff_ID (sort {${$scaff_features{$b}}[1] <=> ${$scaff_features{$a}}[1]} (keys %scaff_features)) {
					
					# Count number of checked scaffolds
					++$L90_genes;
					
					# Add number of genes on this scaffold to checked gene number
					if ($genes_on_scaff{$scaff_ID}) {
						$gene_count += $genes_on_scaff{$scaff_ID};
					}
					else {
						$gene_count += 0;
					}
			
					# Check whether 90% of total gene number is already reached
					last if ($gene_count >= ($non_iso_transcript_nr*0.9));
				}
				
		## Format lengths
			my $formatted_ass_length 		= $en->format_number($assembly_length);
			my $formatted_ass_length_noN	= $en->format_number($assembly_length_noN);
			my $perc_Ns						= sprintf ("%.2f" ,($assembly_Ns/$assembly_length)*100);
	
	
######## INTRON distribution a la Roy & Penny 2007 ####
	
	# Required values
		## Number of introns
		## 3n: 			Fraction of introns with a multiple of 3 bases length			(count/intron #)
		## 3n+1:			Fraction of introns with a multiple of 3 bases length + 1	(count/intron #)
		## 3n+2:			Fraction of introns with a multiple of 3 bases length + 2	(count/intron #)
		## Excess 3n:	Fraction of 3n-introns deviating from the other length types	(3n - ((3n+1) + (3n+2)/2)/intron #)
		## DevExp 3n:	Fraction of 3n-introns deviating from the expected 3n value		(3n - (number of introns/3)/intron #)
	
	# Get values
		## Prepare variables
			my $nr_introns 	= scalar @{$intron_columns[1]};
			my $three_n 	= 0;
			my $three_n1	= 0;
			my $three_n2	= 0;
	
		## Go through intron lengths
			foreach my $length (@{$intron_columns[1]}) {
				# Check with modulo the lengths
				++ $three_n 	if ($length % 3 == 0);
				++ $three_n1 	if ($length % 3 == 1);
				++ $three_n2	if ($length % 3 == 2);
			}
		
		## Calculate values
			my $excess_3n 	= sprintf "%.3f", ( ($three_n - (($three_n1 + $three_n2)/2)) / $nr_introns);
			my $dev_exp_3n	= sprintf "%.3f", ( ($three_n - ($nr_introns/3)) / $nr_introns);
		

######## Overlap check ###########################
my %overlapping_transcripts;

foreach my $trnscrpt (keys %transcript_loci) {
	${$transcript_loci{$trnscrpt}}[0] =~ /^(.*?):(\d*)-(\d*)$/;
	my $loc_id 	= $1;
	my $start 	= $2;
	my $stop 	= $3;
	my $strand	= ${$transcript_loci{$trnscrpt}}[1];
	
	foreach my $pot_overlap (keys %transcript_loci) {
		# next if already checked!
		next if grep (/^$pot_overlap$/, keys %overlapping_transcripts);
		
		${$transcript_loci{$pot_overlap}}[0] =~ /^(.*?):(\d*)-(\d*)$/;
		my $po_loc_id 	= $1;
		my $po_start 	= $2;
		my $po_stop 	= $3;
		my $po_strand	= ${$transcript_loci{$pot_overlap}}[1];
		
		# must be on same scaffold and same strand
		next if $loc_id ne $po_loc_id;
		next if $strand ne $po_strand;
		
		my $keep_flag = 0;
				
		# check for overlap 
		if (my $is_between_start = (sort {$a <=> $b} $po_start, $po_stop, $start)[1] == $start) {
			$keep_flag = 1;
		}
		elsif (my $is_between_stop = (sort {$a <=> $b} $po_start, $po_stop, $stop)[1] == $stop) {
			$keep_flag = 1;
		}
		
		if ($keep_flag and $trnscrpt ne $pot_overlap) {
			push @{$overlapping_transcripts{$trnscrpt}}, $pot_overlap;
		}
	}
}

	
######## PRINT SUMMARY FILE ######################
	if ($print_file[1]) {
		
		print "# Printing summary file...$n";
		
		print $out_summary	"## Tool:${t}$version$n",
							"## Input files:${t}Gff${t}$feature_file$n",
							"##${t}Fasta${t}$fasta_file\n$n";
		
		### Counting types ###
		print "## Counting types...$n";
		
		my $scaff_count = scalar @scaff_IDs;
		print $out_summary "$n# COUNTS$n";
		print $out_summary "$t${t}Total count (present in annotation)${t}Analyzed (excl. isoforms)${t}%$n";
		print $out_summary "Type: Scaffolds${t}${t}$scaff_count$n";
		
		for my $type (@types) {
			my $type_rs = $features->search({type => $type});
			my %counts;
			# Count types GAL style
			while (my $feature = $type_rs->next) {
				my $children = $feature->children;
				while (my $child = $children->next) {
					$counts{$child->type}++;
				}
			}
			next unless %counts;
			# Print counts for child types
			print $out_summary  "Type: $type$n";
			for my $child_type (keys %counts) {
				print $out_summary  "${t}Child Type: $child_type${t}" . $counts{$child_type};
				
				# Print additional info (how many were analyzed?)
				if 		($type eq 'gene' and $child_type eq 'mRNA') {
					printf $out_summary "${t}%.0f${t}%.2f", $non_iso_transcript_nr,($non_iso_transcript_nr/$counts{$child_type}*100);
					print  $out_summary ' %';
				}
				elsif 	($type eq 'mRNA' and $child_type eq 'intron') {
					printf $out_summary "${t}%.0f${t}%.2f", scalar @{$intron_columns[1]}, (scalar @{$intron_columns[1]}/$counts{$child_type}*100);
					print  $out_summary ' %';
				}
				elsif 	($type eq 'mRNA' and $child_type eq 'exon') {
					printf $out_summary "${t}%.0f${t}%.2f", scalar @{$exon_columns[1]}, (scalar @{$exon_columns[1]}/$counts{$child_type}*100);
					print  $out_summary ' %';
				}
				elsif 	($type eq 'mRNA' and $child_type eq 'CDS') {
					printf $out_summary "${t}%.0f${t}%.2f", scalar @{$cds_columns[1]}, (scalar @{$cds_columns[1]}/$counts{$child_type}*100);
					print  $out_summary ' %';
				}
				
				print $out_summary  "$n";
			}
		}
		print $out_summary $n;
		
		print $out_summary "Single CDS gene count${t}$single_cds_genes$n";
		print $out_summary "Single exon gene count${t}$single_exon_genes$n$n";
		
		print $out_summary "Strand-Mix gene count${t}$strand_mix_genes$n$n";
		
		
		# Assembly values
			print $out_summary "# ASSEMBLY\n";
		
			# Length of assembly
			print $out_summary "Assembly size (bp)${t}${t}$formatted_ass_length$n";
			print $out_summary "Assembly size (bp) without Ns${t}${t}$formatted_ass_length_noN$n";
			print $out_summary "Percentage of Ns in assembly${t}${t}$perc_Ns %$n$n";
		
			# Assembly GCs
			my $formatted_GCs 				= $en->format_number($assembly_GCs);
			print $out_summary "Total number of GC bases in assembly${t}${t}${formatted_GCs}$n";
			print $out_summary "Assembly GC content (GC/total length)${t}${t}${assembly_GCratio} %$n";
			print $out_summary "Assembly GC content without ambiguity (GCS/length-NRYKMBDHV)${t}${t}${assembly_GCnoAmratio} %$n$n";
			
			# Assembly CpG o/e
			print $out_summary "CpG o/e for complete assembly${t}${t}${assembly_CpGoe}$n";
			print $out_summary "CpG o/e for assembly without ambiguity${t}${t}${assembly_CpGoe_noN}$n$n";
			
			# N and L values
			print $out_summary "Values of N50, L50, etc.:$n";
			print $out_summary "N50${t}L50${t}N75${t}L75${t}N90${t}L90$n";
			print $out_summary "$N50${t}$L50${t}$N75${t}$L75${t}$N90${t}$L90$n$n";
			
			print $out_summary "L90pcG${t}$L90_genes$n$n";
			
			# Intron distribution a la Roy & Penny 2007
			print $out_summary "# INTRON LENGTH DISTRIBUTION  (see Roy & Penny 2007)$n";
			print $out_summary "Intron count${t}3n${t}3n+1${t}3n+2${t}Excess 3n${t}Deviation from expected 3n$n";
			print $out_summary "$nr_introns${t}$three_n${t}$three_n1${t}$three_n2${t}$excess_3n${t}$dev_exp_3n$n$n";
		
		
		# Scaffold values 
			print $out_summary "# TRANSCRIPT DATA per Scaffold${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'Transcript count per scaffold', 									$scaff_columns[ 5]);
				print_stats($out_summary, 'Added transcript length per scaffold', 							$scaff_columns[ 6]);
				print_stats($out_summary, 'Transcript coverage (added transcript length / scaff length)', 	$scaff_columns[ 7]);
				print_stats($out_summary, 'Transcript density (transcript count / scaff length)', 			$scaff_columns[ 8]);
				
			print $out_summary "$n# CDS DATA per Scaffold${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n" if $overall_cds_presence;
				print_stats($out_summary, 'CDS count per scaffold', 										$scaff_columns[ 9]) if $overall_cds_presence;
				print_stats($out_summary, 'Added CDS length per scaffold', 									$scaff_columns[10]) if $overall_cds_presence;
				print_stats($out_summary, 'CDS count per scaffold', 										$scaff_columns[11]) if $overall_cds_presence;
				print_stats($out_summary, 'CDS density (CDS count / scaff length)', 						$scaff_columns[12]) if $overall_cds_presence;
				
			print $out_summary "$n# EXON DATA per Scaffold${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n" if $overall_exon_presence;
				print_stats($out_summary, 'Exon count per scaffold',										$scaff_columns[13]) if $overall_exon_presence;
				print_stats($out_summary, 'Added exon length per scaffold', 								$scaff_columns[14]) if $overall_exon_presence;
				print_stats($out_summary, 'Exon coverage (added exon length / scaff length)',				$scaff_columns[15]) if $overall_exon_presence;
				print_stats($out_summary, 'Exon density (exon count / scaff length)', 						$scaff_columns[16]) if $overall_exon_presence;
		
			print $out_summary "$n# INTRON DATA per Scaffold${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'Intron count per scaffold', 										$scaff_columns[17]);
				print_stats($out_summary, 'Added intron length per scaffold', 								$scaff_columns[18]);
				print_stats($out_summary, 'Intron coverage (added intron length / scaff length)', 			$scaff_columns[19]);
				print_stats($out_summary, 'Intron density (intron count / scaff length)', 					$scaff_columns[20]);
		
		
		# Lengths
			print  $out_summary "$n# Individual LENGTHs${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation${t}Sum$n";
				print_stats($out_summary, 'Scaffold',				$scaff_columns[ 1]);
				print_stats($out_summary, 'Transcript (genomic)',	$transcript_columns[ 1], 'sum');	
				print_stats($out_summary, 'Protein', 				$transcript_columns[ 6], 'sum');			
				print_stats($out_summary, 'CDS', 					$cds_columns[ 1], 'sum') if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 					$exon_columns[ 1], 'sum') if $overall_exon_presence;
				print_stats($out_summary, 'Intron', 				$intron_columns[ 1], 'sum');
		
			print  $out_summary "$n# Added LENGTHs per transcript${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 	$transcript_columns[ 9]) 	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[20]) 	if $overall_exon_presence;
				print_stats($out_summary, 'Intron', $transcript_columns[31]);
			
			print  $out_summary "$n# Median LENGTHs per transcript$n";
				print $out_summary "${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 	$transcript_columns[10])	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[21]) 	if $overall_exon_presence;
				print_stats($out_summary, 'Intron',	$transcript_columns[32]);
				
			print  $out_summary "$n# Mean LENGTHs per transcript${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 	$transcript_columns[11]) 	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[22]) 	if $overall_exon_presence;
				print_stats($out_summary, 'Intron', $transcript_columns[33]);
		
		# GC contents
			print $out_summary  "$n# Individual GC CONTENT${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'Transcript',						$transcript_columns[2]);
				print_stats($out_summary, 'Transcript (without ambiguity)',	$transcript_columns[3]);	
				print_stats($out_summary, 'CDS', 							$cds_columns[2])		if $overall_cds_presence;
				print_stats($out_summary, 'CDS (without ambiguity)', 		$cds_columns[3])		if $overall_cds_presence;	
				print_stats($out_summary, 'Exon', 							$exon_columns[2])		if $overall_exon_presence;
				print_stats($out_summary, 'Exon (without ambiguity)',		$exon_columns[3])		if $overall_exon_presence;	
				print_stats($out_summary, 'Intron',							$intron_columns[2]);
				print_stats($out_summary, 'Intron (without ambiguity)', 	$intron_columns[3]);
			
			print $out_summary  "$n# Median GC CONTENT per transcript${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 						$transcript_columns[12])	if $overall_cds_presence;
				print_stats($out_summary, 'CDS (without ambiguity)',	$transcript_columns[14])	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 						$transcript_columns[23])	if $overall_exon_presence;
				print_stats($out_summary, 'Exon (without ambiguity)', 	$transcript_columns[25]) 	if $overall_exon_presence;	
				print_stats($out_summary, 'Intron', 					$transcript_columns[34]);
				print_stats($out_summary, 'Intron (without ambiguity)',	$transcript_columns[36]);
				
			print $out_summary  "$n# Mean GC CONTENT per transcript${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 						$transcript_columns[13])	if $overall_cds_presence;
				print_stats($out_summary, 'CDS (without ambiguity)',	$transcript_columns[15])	if $overall_cds_presence;	
				print_stats($out_summary, 'Exon', 						$transcript_columns[24])	if $overall_exon_presence;
				print_stats($out_summary, 'Exon (without ambiguity)',	$transcript_columns[26])	if $overall_exon_presence;	
				print_stats($out_summary, 'Intron', 					$transcript_columns[35]);
				print_stats($out_summary, 'Intron (without ambiguity)',	$transcript_columns[37]);
		
		# CpG o/e
			print $out_summary  "$n# Individual CpG o/e values${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'Transcript',					$transcript_columns[ 4]);
				print_stats($out_summary, 'CDS', 						$cds_columns[4]) 		if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 						$exon_columns[4])		if $overall_exon_presence;
				print_stats($out_summary, 'Intron',						$intron_columns[4]);
				
			print $out_summary  "$n# Median CpG o/e values per transcript${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 	$transcript_columns[16])	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[27]) 	if $overall_exon_presence;
				print_stats($out_summary, 'Intron',	$transcript_columns[38]);
		
		# Counts	
			print  $out_summary "$n# COUNTS per transcript${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'Alternative Spliceforms', 	$transcript_columns[ 7]);
				print_stats($out_summary, 'CDS', 	$transcript_columns[ 8])	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[19])	if $overall_exon_presence;
				print_stats($out_summary, 'Intron',	$transcript_columns[30]);
		
		# Coverages
			print  $out_summary "$n# COVERAGES per transcript (length x / transcript length)${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 	$transcript_columns[17])	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[28])	if $overall_exon_presence;
				print_stats($out_summary, 'Intron',	$transcript_columns[39]); 
		
		# Densities
			print  $out_summary "$n# DENSITIES per transcript (# x / transcript length)${t}Minimum${t}Maximum${t}Mean${t}Median${t}Variance${t}Standard deviation$n";
				print_stats($out_summary, 'CDS', 	$transcript_columns[18])	if $overall_cds_presence;
				print_stats($out_summary, 'Exon', 	$transcript_columns[29])	if $overall_exon_presence;
				print_stats($out_summary, 'Intron',	$transcript_columns[40]); 
		
		# Overlapping transcripts
			if (keys %overlapping_transcripts) {
				my @lines;
				my $overlap_count;
				
				foreach my $overlapped (sort keys %overlapping_transcripts) {
					$overlap_count += scalar @{$overlapping_transcripts{$overlapped}};
					my $line = join ("\t", $overlapped, @{$overlapping_transcripts{$overlapped}}, $n);
					push @lines, $line;
				}
				
				print $out_summary "$n# OVERLAPPING TRANSCRIPTS (same strand)$n";
				print $out_summary "${t}Count${t}Percentage (Count/Gene count)$n";
				print $out_summary "Overlapping transcripts${t}$overlap_count${t}", ($overlap_count/$non_iso_transcript_nr)*100, $n,$n;
				print "${t}Found $overlap_count overlapping transcripts. Check the summary file to identify them.$n";
				print $out_summary "## Transcript ID${t}Overlapping transcript IDs$n";
				print $out_summary @lines;
				
			}
		
		
		# Missing scaffolds
			if ($scaff_miss_flag) { 
				my @lines;
				my $missed_transcripts;
				
				foreach my $missed_scaff (sort keys %missing_scaffs) {
					$missed_transcripts += scalar @{$missing_scaffs{$missed_scaff}};
					my $line = join ("\t", $missed_scaff, @{$missing_scaffs{$missed_scaff}}, $n);
					push @lines, $line;
				}
				
				print $out_summary "$n# MISSING SCAFFOLDS in fasta, skipped transcripts$n";
				print $out_summary "Missing scaffolds${t}", (scalar @lines), $n;
				print $out_summary "Skipped transcripts${t}$missed_transcripts$n$n";
				print $out_summary "## Scaff ID${t}transcript ID$n";
				print $out_summary @lines;
			}
		
		# Print (STDOUT) missing scaffolds
		print "NOTE: ", scalar keys %missing_scaffs, " scaffolds missing in the fasta (also see summary!):$n",
				join (', ', keys %missing_scaffs),$n if $scaff_miss_flag;
			
	}
	
	
######## PRINT SCAFFOLD FILES #######################
	# Determine whether scaffold data is to be printed at all
		my $print_scaff_flag = 0;
		foreach my $file_nr (2 .. 6) {
			if ($print_file[$file_nr] == 1) {
				$print_scaff_flag = 1;
				last;
			}
		}
	
	# If yes, print scaffold data
	if ($print_scaff_flag) {
		print "# Printing scaffold data...$n";
		
		# Prepare header
			my $head_s_general	= "Scaffold length${t}Scaffold length without Ns${t}Scaffold GC content${t}".
									"Scaffold GC content without ambiguity${n}";
			my $head_s_t		= "Transcript count of scaffold${t}Added transcript length of scaffold${t}".
									"Transcript coverage of scaffold${t}Transcript density of scaffold${t}".
									"Transcripts on + : - strand${t}${n}";
			my $head_s_cds		= "CDS count of scaffold${t}Added CDS length of scaffold${t}CDS coverage of scaffold${t}".
									"CDS density of scaffold${n}";
			my $head_s_e		= "Exon count of scaffold${t}Added exon length${t}Exon coverage of scaffold${t}".
									"Exon density of scaffold$n";
			my $head_s_i		= "Intron count of scaffold${t}Added intron length of scaffold${t}".
									"Intron coverage of scaffold${t}Intron density of scaffold$n";
		
		# Print ID referenced files
			
			## Print header to desired scaffold files
				print $OUT_s_general "Scaffold ID${t}", $head_s_general	if $print_file[ 2];
				print $OUT_s_t		 "Scaffold ID${t}", $head_s_t		if $print_file[ 3];
				print $OUT_s_cds 	 "Scaffold ID${t}", $head_s_cds		if $print_file[ 4];
				print $OUT_s_e 		 "Scaffold ID${t}", $head_s_e		if $print_file[ 5];
				print $OUT_s_i 		 "Scaffold ID${t}", $head_s_i		if $print_file[ 6]; 
		
		
			# Print scaffold data to desired files (sorted by scaffold ID)
						# first to last row 
				for my $i (0 .. $#scaff_rows_IDsort) {
					
					### Scaffold general
					print $OUT_s_general join ("\t", ${$scaff_rows_IDsort[$i]}[0], @{$scaff_rows_IDsort[$i]}[1..4])	if $print_file[ 2];
					
					### Scaffold transcripts
						print $OUT_s_t 	join ("\t", ${$scaff_rows_IDsort[$i]}[ 0], @{$scaff_rows_IDsort[$i]}[5..8], 
													${$scaff_rows_IDsort[$i]}[23])	if $print_file[ 3];
						
					### Scaffold CDS
					if ($overall_cds_presence) {
						print $OUT_s_cds join ("\t", ${$scaff_rows_IDsort[$i]}[0], @{$scaff_rows_IDsort[$i]}[9..12])	if $print_file[ 4];
					}
					
					### Scaffold exons
					if ($overall_exon_presence) {
						print $OUT_s_e 	join ("\t", ${$scaff_rows_IDsort[$i]}[0], @{$scaff_rows_IDsort[$i]}[13..16])	if $print_file[ 5];
					}
					
					### Scaffold introns
					if (defined ${$scaff_rows_IDsort[$i]}[17]) { # if i_counts for this transcript
						print $OUT_s_i 	join ("\t", ${$scaff_rows_IDsort[$i]}[0], @{$scaff_rows_IDsort[$i]}[17..20])	if $print_file[ 6];
					}
					
					print $OUT_s_general	"$n"	if $print_file[ 2];
					print $OUT_s_t 			"$n"	if $print_file[ 3];
					print $OUT_s_cds 		"$n"	if $print_file[ 4];
					print $OUT_s_e 			"$n"	if $print_file[ 5];
					print $OUT_s_i			"$n"	if $print_file[ 6];
				}
				
	}
	
	
		
	
############### PRINT TRANSCRIPT FILES ######################
	# Determine whether transcript data is to be printed at all
		my $print_transcript_flag = 0;
		foreach my $file_nr (7 .. 10) {
			if ($print_file[$file_nr] == 1) {
				$print_transcript_flag = 1;
				last;
			}
		}
	
	# If yes, print transcript data
	if ($print_transcript_flag) {
	
		print "# Printing transcript data...$n";	
		
		
		# Prepare headers
			my $head_t_general	= "Transcript length (genomic)${t}Transcript GC content${t}Transcript GC content without ambiguity${t}".
									"Transcript CpG o/e${t}Transcript strand${t}Protein length${t}Count of alternative spliceforms${n}";
			my $head_t_cds		= "CDS count of transcript${t}Added CDS length of transcript${t}".
									"Median CDS length of transcript${t}Mean CDS length of transcript${t}".
									"Median CDS GC content of transcript${t}Mean CDS GC content of transcript${t}".
									"Median CDS GC content without ambiguity of transcript${t}Mean CDS GC content without ambiguity of transcript${t}".
									"Median CDS CpG o/e${t}".
									"CDS coverage of transcript${t}CDS density of transcript${n}";
			my $head_t_e		= "Exon count of transcript${t}Added exon length of transcript${t}".
									"Median exon length of transcript${t}Mean exon length of transcript${t}".
									"Median exon GC content of transcript${t}Mean exon GC content of transcript${t}".
									"Median exon GC content without ambiguity of transcript${t}Mean exon GC content without ambiguity of transcript${t}".
									"Median exon CpG o/e${t}".
									"Exon coverage of transcript${t}Exon density of transcript${n}";
			my $head_t_i		= "Intron count of transcript${t}Added intron length of transcript${t}".
									"Median intron length of transcript${t}Mean intron length of transcript${t}".
									"Median intron GC content of transcript${t}Mean intron GC content of transcript${t}".
									"Median intron GC content without ambiguity of transcript${t}Mean intron GC content without ambiguity of transcript${t}".
									"Median intron CpG o/e${t}".
									"Intron coverage of transcript${t}Intron density of transcript$n";
		
		# Print ID references files
			## Print header to desired transcript files
				print $OUT_t_general	"Transcript ID${t}", $head_t_general	if $print_file[ 7];
				print $OUT_t_cds		"Transcript ID${t}", $head_t_cds		if $overall_cds_presence and $print_file[ 8];
				print $OUT_t_e			"Transcript ID${t}", $head_t_e			if $overall_exon_presence and $print_file[ 9];
				print $OUT_t_i			"Transcript ID${t}", $head_t_i			if $print_file[10];
		
			## Print transcript data to desired files (sorted by transcript ID)
				for my $i (0 .. $#transcript_rows_IDsort) {
				
					### t_general
					print $OUT_t_general join ("\t", ${$transcript_rows_IDsort[$i]}[0], @{$transcript_rows_IDsort[$i]}[1..7])	if $print_file[ 7];
					
					### t_cds
					if ($overall_cds_presence) {
						print $OUT_t_cds join ("\t", ${$transcript_rows_IDsort[$i]}[0], @{$transcript_rows_IDsort[$i]}[8..18])	if $print_file[ 8];
					}
					
					### t_e
					if ($overall_exon_presence) {
						print $OUT_t_e join ("\t", ${$transcript_rows_IDsort[$i]}[0], @{$transcript_rows_IDsort[$i]}[19..29])	if $print_file[ 9];
					} 
					
					### t_i
					if (defined ${$transcript_rows_IDsort[$i]}[32]) { # if i_counts for this transcript
						print $OUT_t_i join ("\t", ${$transcript_rows_IDsort[$i]}[0], @{$transcript_rows_IDsort[$i]}[30..40])	if $print_file[10];
					}
					
					print $OUT_t_general	"$n"	if $print_file[ 7];
					print $OUT_t_cds		"$n"	if $print_file[ 8];
					print $OUT_t_e			"$n"	if $print_file[ 9];
					print $OUT_t_i			"$n"	if $print_file[10];
				}
	
	}
	
	
	
###################### CDS, Exon, Intron files ########################
	
	# Prepare header
		my $head_cds_e_i = "ID${t}Transcript ID${t}Length${t}GC content${t}GC content without ambiguity${t}CpG o/e${t}Strand$n";
	
	# Print CDS file
		if ($print_file[11]) {
			print "# Printing CDS data...$n";	
			print $OUT_cdss $head_cds_e_i;
			foreach my $ID (sort keys %cds_features) {
				print $OUT_cdss join ($t, $ID, @{$cds_features{$ID}}), $n;
			}
		}
		
	# Print exon file
		if ($print_file[12]) {
			print "# Printing exon data...$n";	
			print $OUT_exons $head_cds_e_i;
			foreach my $ID (sort keys %exon_features) {
				print $OUT_exons join ($t, $ID, @{$exon_features{$ID}}), $n;
			}
		}
	
	# Print intron file
		if ($print_file[13]) {
			print "# Printing intron data...$n";	
			print $OUT_introns $head_cds_e_i;
			foreach my $ID (sort{$a <=> $b} keys %intron_features) {
				print $OUT_introns join ($t, $ID, @{$intron_features{$ID}}), $n;
			}
		}
	
	
###################### SUMMARY LINE to BATCH FILE #####################
	
	if ($batch_name) {
		print "# Printing data to batch files...$n";
		
		# Batch GENERAL
		if ($print_file[14]) {
			my $head_batch_general	=	"Species${t}".
										"Assembly size${t}".
										"Assembly size without Ns${t}".
										"Ns in assembly (%)${t}". 
										"GCs in assembly (total count)${t}".
										"Assembly GC content${t}".
										"Assembly GC content without ambiguity${t}".
										"CpG o/e for complete assembly${t}".
										"CpG o/e for assembly without ambiguity${t}".
										"N50${t}L50${t}N75${t}L75${t}N90${t}L90${t}". 
										"L90(genes)${t}". 
										"Scaffold count (total)${t}".
										"Gene count (total)${t}".
										"CDS count (total)${t}".
										"Exon count (total)${t}".
										"Intron count (total)${t}".
										"3n intron count${t}3n+1 intron count${t}".
										"3n+2 intron count${t}Excess 3n intron fraction${t}Deviation from expected 3n intron fraction${n}";
										
			## If file is empty, print header
			if (!-s $OUT_batch_general) {
				print $OUT_batch_general $head_batch_general;
			}
			## Append data line
			print $OUT_batch_general join("\t", $out_name, 
												$formatted_ass_length, 
												$formatted_ass_length_noN, 
												$perc_Ns,
												$assembly_GCs,
												$assembly_GCratio, 
												$assembly_GCnoAmratio, 
												$assembly_CpGoe, $assembly_CpGoe_noN,
												$N50, $L50, $N75, $L75, $N90, $L90,
												$L90_genes,
												scalar @scaff_IDs,
												$non_iso_transcript_nr,
												scalar @{$cds_columns[ 1]},
												scalar @{$exon_columns[ 1]},
												scalar @{$intron_columns[ 1]},
												$three_n, $three_n1, $three_n2, $excess_3n, $dev_exp_3n,
										), $n;
		}
	
		# Batch SCAFFOLD means
		
		if ($print_file[15]) {
			my $head_batch_s_means	=	"Species${t}Mean scaffold length$t".
										"Mean scaffold length -Ns$t".
										"Mean scaffold GC content$t".
										"Mean scaffold GC content without ambiguity$t".
										"Mean transcript count per scaffold$t".
										"Mean added transcript length per scaffold$t".
										"Mean transcript coverage per scaffold$t".
										"Mean transcript density per scaffold$t".
										"Mean CDS count per scaffold$t".
										"Mean added CDS length per scaffold$t".
										"Mean CDS coverage per scaffold$t".
										"Mean CDS density per scaffold$t".
										"Mean exon count per scaffold$t".
										"Mean added exon length per scaffold$t".
										"Mean exon coverage per scaffold$t".
										"Mean exon density per scaffold$t".
										"Mean intron count per scaffold$t".
										"Mean added intron length per scaffold$t".
										"Mean intron coverage per scaffold$t".
										"Mean intron density per scaffold$n";
			
			## If file is empty, print header
				if (!-s $OUT_batch_s_means) {
					print $OUT_batch_s_means $head_batch_s_means;
				}
			## Print data line, start with species name
				print $OUT_batch_s_means $out_name, $t;
			## Continue with values
				print_batch_data($OUT_batch_s_means, 'mean', \@scaff_columns, [1 .. 20]);
			# Finish line with line break
				print $OUT_batch_s_means $n;
		
		}
		
		
		# Batch SCAFFOLD medians
		
		if ($print_file[16]) {
			my $head_batch_s_medians =	"Species${t}Median scaffold length$t".
										"Median scaffold length -Ns$t".
										"Median scaffold GC content$t".
										"Median scaffold GC content without ambiguity$t".
										"Median transcript count per scaffold$t".
										"Median added transcript length per scaffold$t".
										"Median transcript coverage per scaffold$t".
										"Median transcript density per scaffold$t".
										"Median CDS count per scaffold$t".
										"Median added CDS length per scaffold$t".
										"Median CDS coverage per scaffold$t".
										"Median CDS density per scaffold$t".
										"Median exon count per scaffold$t".
										"Median added exon length per scaffold$t".
										"Median exon coverage per scaffold$t".
										"Median exon density per scaffold$t".
										"Median intron count per scaffold$t".
										"Median added intron length per scaffold$t".
										"Median intron coverage per scaffold$t".
										"Median intron density per scaffold$n";
										
			## If file is empty, print header
				if (!-s $OUT_batch_s_medians) {
					print $OUT_batch_s_medians $head_batch_s_medians;
				}		
			## Print data line, start with species name				
				print $OUT_batch_s_medians $out_name, $t;
			## Continue with values
				print_batch_data($OUT_batch_s_medians, 'median', \@scaff_columns, [1 .. 20]);
			## Finish line with line break
				print $OUT_batch_s_medians $n;
		}
		
		# Batch TRANSCRIPT means
		if ($print_file[17]) {
			
			my $head_batch_t_means	=	"Species${t}Mean genomic transcript length$t".
										"Mean transcript GC content$t".
										"Mean transcript GC content without ambiguity$t".
										"Mean transcript CpG o/e$t".
										"Mean protein length$t".
										"Mean count of alternative spliceforms$t".
										"Mean individual CDS length$t".
										"Mean individual CDS GC content$t".
										"Mean individual CDS GC content without ambiguity$t".
										"Mean individual exon length$t".
										"Mean individual exon GC content$t".
										"Mean individual exon GC content without ambiguity$t".
										"Mean individual intron length$t".
										"Mean individual intron GC content$t".
										"Mean individual intron GC content without ambiguity$t".
										"Mean CDS count per transcript$t".
										"Mean added CDS length per transcript$t".
										"Mean median CDS length per transcript$t".
										"Mean median CDS GC content per transcript$t".
										"Mean median CDS GC content without ambiguity per transcript$t".
										"Mean median CDS CpG o/e per transcript$t".
										"Mean CDS coverage per transcript$t".
										"Mean CDS density per transcript$t".
										"Mean exon count per transcript$t".
										"Mean added exon length per transcript$t".
										"Mean median exon length per transcript$t".
										"Mean median exon GC content per transcript$t".
										"Mean median exon GC content without ambiguity per transcript$t".
										"Mean median exon CpG o/e per transcript$t".
										"Mean exon coverage per transcript$t".
										"Mean exon density per transcript$t".
										"Mean intron count per transcript$t".
										"Mean added intron length per transcript$t".
										"Mean median intron length per transcript$t".
										"Mean median intron GC content per transcript$t".
										"Mean median intron GC content without ambiguity per transcript$t".
										"Mean median intron CpG o/e per transcript$t".
										"Mean intron coverage per transcript$t".
										"Mean intron density per transcript$n";
			
			## If file is empty, print header
				if (!-s $OUT_batch_t_means) {
					print $OUT_batch_t_means $head_batch_t_means;
				}
			## Print data line, start with species name
				print $OUT_batch_t_means $out_name, $t;
			## Continue with values
				print_batch_data($OUT_batch_t_means, 'mean', \@transcript_columns, 	[1 .. 4, 6, 7]);
				print $OUT_batch_t_means $t;
				print_batch_data($OUT_batch_t_means, 'mean', \@cds_columns, 		[1 .. 3]);
				print $OUT_batch_t_means $t;
				print_batch_data($OUT_batch_t_means, 'mean', \@exon_columns, 		[1 .. 3]);
				print $OUT_batch_t_means $t;
				print_batch_data($OUT_batch_t_means, 'mean', \@intron_columns, 		[1 .. 3]);
				print $OUT_batch_t_means $t;
				print_batch_data($OUT_batch_t_means, 'mean', \@transcript_columns, 	[8 .. 10, 12, 14, 16 .. 21, 23, 25, 27 .. 32, 34, 36, 38 .. 40]);
			# Finish line with line break
				print $OUT_batch_t_means $n;
			
		}
		
		
		# Batch TRANSCRIPT medians
		
		if ($print_file[18]) {
			
			my $head_batch_t_medians =	"Species${t}Median genomic transcript length$t".
										"Median transcript GC content$t".
										"Median transcript GC content without ambiguity$t".
										"Median transcript CpG o/e$t".
										"Median protein length$t".
										"Median count of alternative spliceforms$t".
										"Median individual CDS length$t".
										"Median individual CDS GC content$t".
										"Median individual CDS GC content without ambiguity$t".
										"Median individual exon length$t".
										"Median individual exon GC content$t".
										"Median individual exon GC content without ambiguity$t".
										"Median individual intron length$t".
										"Median individual intron GC content$t".
										"Median individual intron GC content without ambiguity$t".
										"Median CDS count per transcript$t".
										"Median added CDS length per transcript$t".
										"Median median CDS length per transcript$t".
										"Median median CDS GC content per transcript$t".
										"Median median CDS GC content without ambiguity per transcript$t".
										"Median median CDS CpG o/e per transcript$t".
										"Median CDS coverage per transcript$t".
										"Median CDS density per transcript$t".
										"Median exon count per transcript$t".
										"Median added exon length per transcript$t".
										"Median median exon length per transcript$t".
										"Median median exon GC content per transcript$t".
										"Median median exon GC content without ambiguity per transcript$t".
										"Median median exon CpG o/e per transcript$t".
										"Median exon coverage per transcript$t".
										"Median exon density per transcript$t".
										"Median intron count per transcript$t".
										"Median added intron length per transcript$t".
										"Median median intron length per transcript$t".
										"Median median intron GC content per transcript$t".
										"Median median intron GC content without ambiguity per transcript$t".
										"Median median intron CpG o/e per transcript$t".
										"Median intron coverage per transcript$t".
										"Median intron density per transcript$n";
			
			## If file is empty, print header
				if (!-s $OUT_batch_t_medians) {
					print $OUT_batch_t_medians $head_batch_t_medians;
				}
			## Print data line, start with species name
				print $OUT_batch_t_medians $out_name, $t;
			## Continue with values
				print_batch_data($OUT_batch_t_medians, 'median', \@transcript_columns, 	[1 .. 4, 6, 7]);
				print $OUT_batch_t_medians $t;
				print_batch_data($OUT_batch_t_medians, 'median', \@cds_columns, 		[1 .. 3]);
				print $OUT_batch_t_medians $t;
				print_batch_data($OUT_batch_t_medians, 'median', \@exon_columns, 		[1 .. 3]);
				print $OUT_batch_t_medians $t;
				print_batch_data($OUT_batch_t_medians, 'median', \@intron_columns, 		[1 .. 3]);
				print $OUT_batch_t_medians $t;
				print_batch_data($OUT_batch_t_medians, 'median', \@transcript_columns, 	[8 .. 10, 12, 14, 16 .. 21, 23, 25, 27 .. 32, 34, 36, 38 .. 40]);
			# Finish line with line break
				print $OUT_batch_t_medians $n;
		
		}
	
		
		# COMPONENT SIZES
		
		if ($print_file[19]) {
			## If file is empty, print header	
				if (!-s $out_compsize) {
					print $out_compsize "Species${t}",
									"Assembly size (Mb)${t}",
									"Coding amount (Mb)${t}",
									"Intron amount (Mb)${t}",
									"Transposable elements (Mb)${t}",
									"Other repeats (Mb)${t}",
									"Rest size (Mb)${t}${t}",
									"Coding amount (%)${t}",
									"Intron amount (%)${t}",
									"Transposable elements (%)${t}",
									"Other repeats (%)${t}",
									"Rest size (%)$n";	
				}
			
			## Calculate values, format as MB / %
				my $assembly_size_Mb = $assembly_length/1000000;
				my $CDS_size_Mb 	 = (sum_def(@{$cds_columns[ 1]}))/1000000;
				my $i_size_Mb		 = (sum_def(@{$intron_columns[ 1]}))/1000000;
				my $perc_CDS		 = sprintf ("%.3f", ($CDS_size_Mb*100)/$assembly_size_Mb);
				my $perc_i			 = sprintf ("%.3f", ($i_size_Mb*100)/$assembly_size_Mb);
			
			## Prepare line for printing
				my $compsize_line = join("\t", $out_name, $assembly_size_Mb, $CDS_size_Mb, $i_size_Mb, '', '', '', '', $perc_CDS,$perc_i, '', '', ''
									);
			## Print
			print $out_compsize $compsize_line, $n;
		}
		
		# BASH COMMAND FILE
		
		if ($print_file[20]) {
			my $command = "
echo '$out_name'
python BUSCO_v1.1.py \\
-in ${cwd}/COGNATE_${out_name}/${out_name}_00-analyzed_transcripts.fa \\
-o $out_name \\
-l A -m ogs -f \\
1>> busco_LOG.log 2>&1 
";
		
			print $out_bash $command, $n;
		}
	}
	
	
	
	
	
	
###################### FINISH ##########################################
	close $out_fa				if $print_file[ 0];
	close $out_summary			if $print_file[ 1];
	close $OUT_s_general		if $print_file[ 2];
	close $OUT_s_t				if $print_file[ 3];
	close $OUT_s_cds			if $print_file[ 4];
	close $OUT_s_e				if $print_file[ 5];
	close $OUT_s_i				if $print_file[ 6];
	close $OUT_t_general		if $print_file[ 7];
	close $OUT_t_cds			if $print_file[ 8];
	close $OUT_t_e				if $print_file[ 9];
	close $OUT_t_i				if $print_file[10];
	close $OUT_cdss				if $print_file[11];
	close $OUT_exons			if $print_file[12];
	close $OUT_introns			if $print_file[13];
	close $OUT_batch_general	if $batch_name and $print_file[14];
	close $OUT_batch_s_means	if $batch_name and $print_file[15];
	close $OUT_batch_s_medians	if $batch_name and $print_file[16];
	close $OUT_batch_t_means	if $batch_name and $print_file[17];
	close $OUT_batch_t_medians	if $batch_name and $print_file[18];
	close $out_compsize			if $batch_name and $print_file[19];
	close $out_bash				if $batch_name and $print_file[20];
	  

	# Local time stamp
	my ($second, $minute, $hour, $day, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year 		= 1900 + $yearOffset;
	$month 			= $month +1;
	my $theTime 	= sprintf "$year-$month-$day %02d:%02d:%02d", $hour, $minute, $second;
	
	# Runtime
	$now 		= time - $now;
	my $runTime = sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);

	print "$n## DONE with $out_name! ($theTime, run time: $runTime)$n$n";
	
	
	
	
	# Return to working dir (outside output dir)
	chdir $cwd;

}





########################################################################
#----------------------------------------------------------------------#
#                     APPENDIX                                         #
#----------------------------------------------------------------------#
################## SUBS ################################################
=head1 SUBROUTINES

=cut

###############################################
=head2 set_print_files

 Usage   : @print_file = set_print_files(\@print, \@dont_print);
 Function: Takes user input which files tp print and which not, stores in array for late ruse by tool (file switches).
 Returns : @print_file (0 or 1 for each file)
 Args    : References of print and dont_print options
 
=cut

sub set_print_files {

	my $print_REF 		= shift @_;
	my $dont_print_REF	= shift @_;
	

	my @print 		= @{$print_REF};
	my @dont_print 	= @{$dont_print_REF};
	
	## Default: print all files
	my	@print_file =  (	1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
							1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	
	# Do not print files not in @print if @print is used
	if (@print) {
		## Prepare list
		my @print = split (/,\s?/, join (',', @print));
		
		## Switch off all files
		@print_file =  (	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		
		## Go through list of desired files
		foreach my $i (@print) {
			### Check whether indices are between 0 and 20
			die "--print LIST may contain only values between 0 and 20! Found value ${i}.$n$usage" if $i < 0 or $i > 20;
			
			### Switch on desired files
			$print_file[$i] = 1;
		}
	}
	
	# Do not print files of @dont_print if @dont_print is used
	if (@dont_print) {
		## Prepare list
		my @dont_print = split (/,\s?/, join (',', @dont_print));
		
		## Go through list of undesired files
		foreach my $i (@dont_print) {
			### Check whether indices are between 0 and 20
			die "--dont_print LIST may contain only values between 0 and 20! Found value ${i}.$n$usage" if $i < 0 or $i > 20;
		
			### Switch off undesired files
			$print_file[$i] = 0;
		}
	}
	
	return @print_file;
	
}

###############################################
=head2 get_input_data

 Usage   : %input = get_input_data($input_file);
 Function: Checks and slurps COGNATE.input file, assigns input linewise to
			storage hash %input (name => $feature_file, $fasta_file).
 Returns : %input
 Args    : $input_file
 
=cut

sub get_input_data {
	
	my $input = shift @_;
	
	# Check whether input file is usable
	die "FATAL: Input file does not exist ($input): $!$n" if !-e $input;
	die "FATAL: Input file must be a COGNATE.input file.$n$synopsis" if $input !~ /COGNATE.input$/;
	
	# Prepare variable
	## %input = species => (.gff3, .fa)
	my %input;
	
	# Slurp input file linewise
	my @input_lines = slurp_file($input);

	# Go through lines and extract pairs (.fa and -gff3)
	foreach my $line (@input_lines) {
		
		## Skip comment lines
		next if $line =~ /^#/;
		
		## Split line, get items
		my @items = split ("\t", $line);
		die "FATAL: Input file format not correct (has to contain three items, tab separated). Reformat and submit again!",
			"${n}Items:$n - ", join ("$n - ", @items), "$n$synopsis" if scalar @items != 3;
		
		my ($feature_file, $fasta_file) = assign_files(($items[1], $items[2]));
		
		## Put name as key and files as value into hash
		if (!exists $input{$items[0]}) {
			@input{$items[0]} = [$feature_file, $fasta_file];
		}
		else {
			die "FATAL: Input file may not contain lines with the same name! (Offending: '$items[0]')${n}Reformat and submit again!$n$synopsis";
		}
	}

	return %input;
}

###############################################
=head2 set_outname

 Usage   : $out_name = set_outname($out_name);
 Function: Checks whether given $out_name is defined and if so, whether it is white space-free. 
			If $out_name is not defined (-n was not used), a default name (genome_#) will be set
			and automatically incremented for the next call.
 Returns : $out_name
 Args    : $out_name
 
=cut

sub set_outname {
	my $out_name = shift @_;
	
	# If a name is defined by the user
	if(defined($out_name)){
		
		# Empty string shall be replaced by default string
		if ($out_name eq '') {
			$out_name = undef;
		}
		
		# Replace white spaces with underlines
		elsif ($out_name =~ /[\s]/){
			warn "NOTE: Name contains invalid white space characters. Will replace white spaces with underline character '_'$n";
			$out_name =~ s/\s/\_/g;
		}
	}
	# If no name was set
	if(!defined($out_name)){
		# Number as suffix for default name
		my $no = 1;
		# Flag to keep it
		my $take_it = 0;
		
		while ($take_it == 0) {
			# Generate a default name
			$out_name = "genome_$no";
			# Check whether it already exists
			if(!-e "COGNATE_$out_name"){
				# If not, take this name
				$take_it = 1;
				last;
			}else{
				# Else, try the next one
				$no++;
			}
		}
		warn "NOTE: No name was set. Will use $out_name.\n";
	}
	return $out_name;
}

###############################################
=head2 assign_files

 Usage   : ($feature_file, $fasta_file) = assign_files(@file_list);
 Function: Checks which file is a gff and which is a fasta (throws error if given files do not
			include one gff and one fasta) and assigns them to the respective variables.
			Checks whether paths are absolute (required by GAL) and whether files exist.
 Returns : path and name of gff and fasta files (input)
 Args    : @file_list
 
=cut

sub assign_files {
	my @inputfiles = @_;
	
	my $fa_suffix = '\.fa|\.fn|\.fasta|\.fas|\.fna';
	
	# Check if both files have the same extension
	if ($inputfiles[0] =~ /\.gff.*$/ and $inputfiles[1] =~ /\.gff.?$/) {
		die "FAIL: Two GFF files, it may only be one!
		File names: $inputfiles[0]
				$inputfiles[1]
		$synopsis";
	}
	elsif ($inputfiles[0] =~ /($fa_suffix)$/ and $inputfiles[1] =~ /($fa_suffix)$/) {
		die "FAIL: Two fasta files, it may only be one!
		File names: $inputfiles[0]
				$inputfiles[1]
		$synopsis";
	}
	
	my ($feature_file, $fasta_file);

	# Check first file
		## If first file is a gff...
		if ($inputfiles[0] =~ /\.gff.?$/){
			$feature_file = $inputfiles[0];
		}
		## or a fasta...
		elsif ($inputfiles[0] =~ /($fa_suffix)$/){
			$fasta_file = $inputfiles[0];
		}
		# if neither gff nor fasta...
		else {
			die "FAIL: Input file '$inputfiles[0]' is neither gff(3) nor fasta, judging from the suffix.$n$synopsis";
		}
	# Check second file
		if ($inputfiles[1] =~ /\.gff.?$/){
			$feature_file = $inputfiles[1];
		}
		elsif ($inputfiles[1] =~ /($fa_suffix)$/){
			$fasta_file = $inputfiles[1];
		}
		else {
			die "FAIL: Input file '$inputfiles[1]' is neither gff(3) nor fasta, judging from the suffix.$n$synopsis";
		}
	
	# Check file existence and extensions
	die "FAIL: Paths to files must be absolute! ($feature_file, $fasta_file)$n" if $feature_file =~ /^.\// or $fasta_file =~ /^.\//;
	die "FAIL: Gff file '$feature_file' not found!$n$synopsis" if !-e $feature_file;
	die "FAIL: Fasta file '$fasta_file' not found!$n$synopsis" if !-e $fasta_file;
	
	return $feature_file, $fasta_file;
}

###############################################
=head2 goto_workingdir

 Usage   : goto_workingdir($workingdir);
 Function: Checks whether given path (working dir) is an existing directory. If yes, 
			changes to it.
 Returns : -
 Args    : /PATH/
 
=cut

sub goto_workingdir {
	
	my $workingdir = shift @_;
	
	# Check if exists
	die "FAIL: Working dir $workingdir does not exist.$n$usage" if !-d $workingdir;
	
	# Go to dir
	chdir $workingdir;
	
	return 1;
}

###############################################
=head2 remove_dir

 Usage   : remove_dir($out_name, $def_overwrite);
 Function: Checks whether output directory for given $out_name already exists and, depending on 
			the user's choice, deletes this directory or not.
 Returns : -
 Args    : $out_name, $def_overwrite
 
=cut

sub remove_dir {
	my $out_name 		= shift @_;
	my $def_overwrite 	= shift @_;
	
	# If output directory already exists, remove it?
	if (-e "COGNATE_$out_name") {
		
		# User did not choose default overwrite
		if (!$def_overwrite) {
			
			# Ask whether it shall be overwritten
			print "# Dir exists. overwrite? [y/n]$n";
			my $overwrite = <STDIN>;
			
			# Evaluate answer
			die "STOP: You quit. Start again?$n" if ($overwrite =~ /quit|exit/i);
			die "STOP: Answer incomprehensible. Start again.$n" if ($overwrite !~ /y|yes|n|no/i);
			die "STOP: Do not overwrite -> change --name STRING and start again!$n" if ($overwrite =~ /n|no/i);
			remove_tree("COGNATE_$out_name") if ($overwrite =~ /y|yes/i);
		}
		
		# User says 'overwrite!' via option
		elsif ($def_overwrite) {
			print "# Dir for $out_name exists but is overwritten automatically$n";
			remove_tree("COGNATE_$out_name");
		}			
	}
	
	return 1;
}

###############################################
=head2 goto_outputdir

 Usage   : goto_outputdir($out_name);
 Function: Creates directory (COGNATE_$out_name) and changes to it.
 Returns : -
 Args    : $out_name
 
=cut

sub goto_outputdir {
	# We ensured already, output directory does not exist
	my $out_name = shift @_;
	
	# Make output dir
	mkdir "COGNATE_$out_name" or die "Error: Cannot make dir COGNATE_$out_name: $!$n" ;
	print "# Made dir for $out_name $n";
	
	# Go to output dir
	chdir "COGNATE_$out_name" or die "Error: Cannot change to dir COGNATE_$out_name: $!$n" ;

	return 1;
}

###############################################
=head2 open_FH

 Usage   : open_FH($out_name, $file_name, $wa, $suffix);
 Function: Open file handle for file '$outname_$file_name.$suffix' to write ($wa = '>') or append ('>>').
 Returns : $filehandle
 Args    : $out_name, file_name, write or append operator, file suffix
 
=cut

sub open_FH {
	# Get specifics what file shall be opened
	my $out_name 	= shift @_;
	my $file_name 	= shift @_; # e.g. "summary"
	my $wa			= shift @_; # write or append? -> ">" or ">>"
	my $suffix		= shift @_; # "tsv" or "txt" or whatever
	
	# Create file name
	my $file = $out_name.'_'.$file_name.'.'.$suffix;
	
	# Open file
	open (my $FH, $wa, $file) or die "ERROR: Cannot open $wa file $file: $!$n";
	
	return $FH;
}

###############################################
=head2 fasta_validity

 Usage   : \%headers = fasta_validity($fasta_file);
 Function: Checks whether given fasta is DNA, containing only IUPAC 
			characters. Uses slurp_fasta. 
 Returns : \%headers (header => 0; shortened headers, only first part of header until whitespace)
 Args    : /PATH/file.fa
 
=cut

sub fasta_validity {
	
	my $fasta = shift @_;
	
	print "# Checking validity of fasta file: $fasta\n";
	
	my %headers;
	
	# Slurp fasta (get sequences and headers)
	my ($seqs_REF, $heads_REF) = slurp_fasta($fasta);
	
	# Go through all sequences (referred to by headers)
	foreach my $header (@$heads_REF) {
		
		# Get seq_ID (first part of header in NCBI genome files)
		if ($header =~ m/(^[\w\.\-]*)\s.*/) {
			my $short_head = $1;
			die "CHECK FATALITY: Header $short_head appears more than once in input fasta $fasta!$n",
				"Check your files!$n" if exists $headers{$short_head};
			$headers{$short_head} = 0;
		}
		else { #no NCBI header
			die "CHECK FATALITY: Header $header appears more than once in input fasta $fasta!$n",
				"Check your files!$n" if exists $headers{$header};
			$headers{$header} = 0;
		}
		
		# Test for alphabet of sequence
		my $seq = $seqs_REF->{$header};
		
			## Diagnose protein
				### valid symbols in protein fasta: A C D E F G H I J K L M N O P Q R S T U V W Y X Z 
				###      which are not in DNA: E F J L P Q Z
			my $count_aa = $seq =~ tr/EFJLPQZefjlpqz//;
			die "CHECK FATALITY: Your fasta $fasta is not valid (sequence '$header' contains diagnostic protein characters, one or more of these: EFJLPQZ)$n" if ($count_aa >= 1);
		
			## Diagnose RNA
				### valid IUPAC symbols in RNA: A B C D G I N O R S T U Y
				###       which are not in DNA: I O
			my $count_rna = $seq =~ tr/IOio//;
			die "CHECK FATALITY: Your fasta $fasta is not valid (sequence '$header' contains diagnostic RNA characters, one or more of these: IO)$n" if ($count_rna >= 1);
			
			## Ensure DNA
				### If neither protein nor RNA, check for only IUPAC characters
				### valid IUPAC symbols in nucleotide DNA fasta: A B C D G H K M N R S T U V W X Y
			my $count_nonIUPAC = $seq =~ tr/ABCDGHKMNRSTUVWXYabcdghkmnrstuvwxy//c;
			### Print offending characters
			$seq =~ tr/ABCDGHKMNRSTUVWXYabcdghkmnrstuvwxy//d if ($count_nonIUPAC >= 1);
			die "CHECK FATALITY: Your fasta $fasta is not valid (sequence '$header' contains non-IUPAC characters: $seq)$n" if ($count_nonIUPAC >= 1);
		
	}
	
	print "# Fasta ok.$n";
	
	return \%headers;
}

###############################################
=head2 gff_validity

 Usage   : gff_validity($gff_file);
 Function: Checks whether given gff is gff3 format (only very basic check, for full validation use your
			favorite gff3_validator (e.g. from GAL), checks whether headers/region IDs of gff and fasta
			match. Uses slurp_file. 
 Returns : -
 Args    : /PATH/file.gff
 
=cut

sub gff_validity {
	my ($gff_file, $headers_REF) = @_;
	
	print "# Checking validity of gff file: $gff_file\n";
	
	# Prepare variables
	my $warn_count = 0;
	my $line_count = 0;
	my %absent_IDs;
	
	# Slurp file (line by line)
	my @lines = slurp_file($gff_file);
	
	# Go through lines
	foreach my $line (@lines) {
		# Check for GFF3 pragma /first line should contain '#gff3')
			## if not present, check externally, then turn internal check off (option -f)
		if ($line_count == 0 and $line !~ /^\#\#gff-version\s*3$/) {
			warn "CHECK WARNING: Gff version pragma is missing! To be sure your gff file is 
			a valid gff3 file, test it using a gff3_validator (e.g. from the
			GAL library)!$n";
			++ $warn_count;
		}
		# Raise line count, so we know we have passed the first line
		++$line_count;
		
		# Ignore further comments
		next if $line =~ /^\#\#/;
		
		# Get line elements (should be tab-separated)
		my @line_elements = split ("\t", $line) or die "CHECK FATALITY: Gff columns not tab separated!Check file!$n";
		my $region_ID = $line_elements[0];
		
		# FORWARD CHECK: Fasta headers present as region IDs in the GFF?
			## If header not exists in fasta, remember missing header
			if (!exists ${$headers_REF}{$region_ID}) {
				++$absent_IDs{$region_ID};
			}
			## If header exists, take note (required for REVERSE CHECK)
			else {++${$headers_REF}{$region_ID};};
		
	}
	
	# Print missing region IDs (Result of FORWARD CHECK)
	my $absent_ID_string = join (" ", keys %absent_IDs);
	if (scalar keys %absent_IDs > 0) {
		warn "CHECK WARNING: Mismatch between GFF and Fasta!${n}The following GFF region ID were not present in the
		fasta:${n}$absent_ID_string$n";
		++ $warn_count;
	}
	
	# REVERSE CHECK: Which headers present in the fasta do not appear as region ID in the GFF?
	my %absent_headers;	
	foreach my $header (keys %{$headers_REF}) {
		++$absent_headers{$header} if ${$headers_REF}{$header} == 0; 
	}
	
	# Print missing headers (Result of REVERSE CHECK)
	my $absent_heads_string = join (" >", keys %absent_headers);
	if (scalar keys %absent_headers > 0) {
		warn "CHECK WARNING: Mismatch between fasta and GFF!${n} The following headers present in the fasta were not present in the
		Gff:${n}$absent_heads_string$n" ;
		++$warn_count;
	}
		
	die "CHECK FAILED: Please check your files or turn off the file check (-f) if you think,
		noted problems are tolerable.$n" if $warn_count > 1;
		
	print "# GFF ok.$n";
		
	return 1;
}

###############################################
=head2 slurp_fasta

 Usage   : $sequence_ofREF, $headersREF = slurp_fasta($fastafile);
 Function: Slurps a given fasta file (uses slurp_file), identifies headers
			and sequences (also across line breaks, then seq is concatenated),
			puts them in a hash and stores headers separately. 
 Returns : \%sequence_of, \@headers
 Args    : /PATH/file.fa
 
=cut

sub slurp_fasta {

	my $file = shift @_; 
	
	# Slurp file (line by line)
	my @lines = slurp_file($file);

	# Prepare variables
	my $header;
	my @headers;
	my %sequence_of;
	
	# Process fasta lines
	foreach my $line ( @lines ) {
		if ($line =~ /^>(.*)$/) {
			$header = $1;
			die "Found two sequences with the same header: $header" if $sequence_of{$header};
			push @headers, $header;
		}
		else {
			# Remove white space
			$line =~ s/\s*//g;
			# Add sequence
			$sequence_of{$header} .= $line;
		}
	}
	
	# Return fasta file content as hash reference
	return \%sequence_of, \@headers;
}

###############################################
=head2 slurp_file

 Usage   : @file_lines = slurp_file($file);
 Function: Returns array of lines in a given (incl PATH) file
 Returns : @file_lines
 Args    : /PATH/file
 
=cut

sub slurp_file {
	my $file = shift @_; 
	
	# Open filehandle
	open my $FILE, '<', $file or die "Couldn't open \"$file\": $!";
	
	# Read in entire file content at once 
	my $file_content = do { local $/; <$FILE> };
	
	# Close filehandle
	close $FILE or die "Could't close \"$file\": $!";
	
	# Convert to unix-style line breaks
	$file_content =~ s/\r\n?/\n/g;
	
	# Split file content on line breaks
	my @file_lines = split ( '\n', $file_content );

	return @file_lines;
}

###############################################
=head2 assembly_NL

 Usage   : ($N, $L) = assembly_NL($lengths_REF, $goal);
 Function: Calculates N and L values for the given list of lenghts. Lengths must be sorted from
			high to low values. The N50/75/90 is the length for which the collection of all 
			scaffolds of that length or longer covers at least 50/75/90% of the total assembly
			size. L50/75/90 is the number of scaffolds equal to or longer than N50/75/90.
 Returns : $N, $L
 Args    : Reference to sorted lengths array, size of goal (e.g., assembly size*0.5)
 
=cut

sub assembly_NL {
	# Get array reference of sorted scaffold lengths
	my $lengths_REF 	= shift @_;
	# and the desired size value
	my $goal			= shift @_;
	
	my $N = 0;
	my $L = 0;
	my $sumN = 0;
	
	foreach my $length (@{$lengths_REF}) {
		# Remember last added length (if goal is reached, this is N)
		$N = $length;
	
		# Add up lengths, beginning with the largest
		$sumN += $length;
		
		# Count number of lengths added
		++ $L;
		
		# If desired value is reached, end and return values
		last if ($N >= $goal);
	}
	
	return $N, $L;
}

###############################################
=head2 get_or_make_ID

 Usage   : ($feature_id, $id_suffix) = get_or_make_ID($transcript_id, $feature, $name, $id_suffix);
 Function: Retrieves feature (CDS or exon) ID from GFF or, if not given there, assigns an ID
			using the given suffix (raised afterwards, so every features gets a unique ID: 
			transcriptID_feat_suffix)
 Returns : $feature_id and suffix
 Args    : transcript ID, feature (CDS or exon), name of the feature (CDS or Exon), id_suffix
 
=cut
sub get_or_make_ID {
			
	my $transcript_id	= shift @_;
	my $feature 		= shift @_;
	my $name			= shift @_;
	my $id_suffix 		= shift @_;
	
	my $feature_id;
	
	# If GAL is able to find an ID (in the GFF)
	if ($feature->feature_id) {
		# assign it
		$feature_id = $feature->feature_id;
	}
	# If there is no ID
	else {
		# create one
		$feature_id = $transcript_id.'_'.$name.'_'.$id_suffix;
		# raise suffix for the next feature without id
		++$id_suffix;
	}
	
	# Return ID and suffix 
	return $feature_id, $id_suffix;
}

###############################################
=head2 get_lengths_and_GCs

 Usage   : ($feature_length, $feature_gc, $feature_gc_noam) = get_lengths_and_GCs($feature);
 Function: Retrieves length and GC content (with and without ambiguity) for a given feature.
 Returns : Length, GC content in %, GC content without ambiguity in %
 Args    : Feature (CDS, exon, or intron)
 
=cut		
sub get_lengths_and_GCs {
	
	my $feature = shift @_;
	
	# Get length
	my $feature_length = $feature->length;
	
	my ($feature_gc, $feature_gc_noam) = 'NA';
	
	# Get GC content
	$feature_gc 		= sprintf ( "%.3f", ($feature->gc_content)*100) 		if $feature->gc_content 		!~ /NA/; 
	$feature_gc_noam	= sprintf ( "%.3f", ($feature->gc_content_noAm)*100) 	if $feature->gc_content_noAm 	!~ /NA/;
	
	# Return values
	return $feature_length, $feature_gc, $feature_gc_noam;
}

###############################################
=head2 sum_def

 Usage   : $sum = sum_def(@array);
 Function: Sums up values of array, removes undef values first. Uses
			List::Util::sum.
 Returns : $sum
 Args    : @array of integers
 
=cut

sub sum_def {
	my @array = @_;

	my $vector	= vector(@array);
	my @values	= grep {defined $_} @array;
	@values		= grep {$_ !~ /^\s$/} @values;
	
	my $sum = sum0(@values);
	
	return $sum;
}

###############################################
=head2 my_mean

 Usage   : $mean = my_mean(\@array);
 Function: Gives mean of values of array, removes undef values first. Uses
			List::Util::mean.
 Returns : $mean
 Args    : @array of integers
 
=cut
	sub my_mean {
		
		my $values_REF = shift @_;
		
		# Take defined values
		my @array 				= grep {defined $_} @{$values_REF};		
		my @array_clean 		= grep {$_ !~ /^NA$/} @array;
		my @array_double_clean	= grep {$_ !~ /^$/} @array_clean;
		
		my $mean = 'NA';
		$mean = mean(@array_double_clean) if scalar @array_double_clean >= 1;
		
		return $mean;		
	}
	

###############################################
=head2 my_median

 Usage   : $median = my_median(@array);
 Function: Gives median of values of array, removes undef values first. Uses
			List::Util::median.
 Returns : $median
 Args    : @array of integers
 
=cut
	sub my_median {
		
		my $values_REF = shift @_;
		
		# Take defined values
		my @array 				= grep {defined $_} @{$values_REF};		
		my @array_clean 		= grep {$_ !~ /^NA$/} @array;
		my @array_double_clean	= grep {$_ !~ /^$/} @array_clean;
		
		my $median = 'NA';
		$median = median(@array_double_clean) if scalar @array_double_clean >= 1;
		
		return $median;		
	}

######################################
=head2 print_stats

 Usage   : print_stats($out_FH, 'Intron', \@array);
 Function: Prints mean, median, variance, stddev, min, and max of a list 
			of values for a given element in tab-style. Uses Statistics::Basic.
 Returns : nothing
 Args    : $out_FH, ELEMENT_NAME, $array_REF

=cut

sub print_stats {
	my ($FH, $type, $values_REF, $sum_flag) = @_;
	
	# Take defined values
	my @array 				= grep {defined $_} @{$values_REF};		
	my @array_clean 		= grep {$_ !~ /^NA$/} @array;
	my @array_double_clean	= grep {$_ !~ /^$/} @array_clean;
	
	# If there are values present
	if (scalar @array_double_clean) {	

		# Calculate statistics
		my $vector		= vector(@array_double_clean);
		
		my $mean 		= sprintf ("%.3f", mean($vector));
		my $variance	= sprintf ("%.3f", variance($vector));
		my $median 		= sprintf ("%.3f", median($vector));
		my $stddev 		= sprintf ("%.3f", stddev($vector));
		my $minimum		= sprintf ("%.3f", min(@array_double_clean));
		my $maximum 	= sprintf ("%.3f", max(@array_double_clean));
		
		my $sum = sum_def(@array_double_clean) if $sum_flag;
		
	
		# Print to designated file
		print $FH  join ("\t", $type, $minimum, $maximum, $mean, $median, $variance, $stddev);
		print $FH $t, $sum if $sum_flag;
		print $FH  "$n";
	}
	
	# No values present
	else {
		# Print this info
		print $FH "$type: empty array...$n$n";
	}
	
	return 1;
	
}

#######################
=head2 transpose

 Usage   : @array_of_arrays_transposed = transpose($array_of_arrays_REF);
 Function: Transposes an array of array, i.e., rows are now columns and vice versa.
 Returns : Transposed AoA
 Args    : $array_of_arrays_REF

=cut
sub transpose{
	my $rows_REF = shift @_;
	
	my @transposed;
	
	# Go through each row
	for my $row (@{$rows_REF}) {
		# and every item within
		for my $column (0 .. $#{$row}) {
			# Put item to transposed position in new array of arrays
			push(@{$transposed[$column]}, $row->[$column]);
		}
	}
	
	return @transposed;
}

#######################
=head2 transpose_sort_transpose

 Usage   : @array_of_arrays_transposed = transpose_sort_transpose($array_of_arrays_REF);
 Function: Transposes an array of array, i.e., rows are now columns and vice versa. Sorts
			columns decreasingly, retransposes AoA.
			EXPERIMENTAL...
 Returns : Transposed AoA
 Args    : $array_of_arrays_REF

=cut

# TODO: transpose_sort_transpose needed?

sub transpose_sort_transpose{
	my $rows_REF = shift @_;
	
	my @transposed_AoA;
	my @sorted_AoA;
	
	#print Dumper $rows_REF;
	
	# Transpose input array of arrays
	# Go through each row
	for my $row (@{$rows_REF}) {
		# and every item within
		for my $column (0 .. $#{$row}) {
			# Put item to transposed position in new array of arrays
			push(@{$transposed_AoA[$column]}, $row->[$column]);
		}
	}
	
	# TODO: try foreach und  sort {$b <=> $a or $b cmp $a} @{$row}
	# TODO: try for (my $i = 1; $i < scalar @$list; $i++) { $list->[$i] = [ sort { $a <=> $b } @{$list->[$i]} ] }
	
	# Sort rows of the transposed array
	for my $row (@transposed_AoA) {
		#print @{$row}, $n;
		#my @sorted;
		#if (@{$row}[0] =~ /^\w\w/) {
		#	print 'header: ', @{$row}, $n;
		#	@sorted = sort {$b cmp $a} @{$row};
		#	$row = \@sorted;
		#	next;
		#}
		#@sorted = sort {$b <=> $a}@{$row};
		#$row = \@sorted;
		
		my @sorted = sort @{$row};
		my @rev_sort = reverse @sorted;
		
		$row = \@rev_sort;
		
	}
	
	print "Sorted $n";
	print Dumper \@transposed_AoA; exit;
	
	# Transpose sorted AoA back
	# Go through each row
	for my $row (@transposed_AoA) {
		# and every item within
		for my $column (0 .. $#{$row}) {
			# Put item to transposed position in new array of arrays
			push(@{$sorted_AoA[$column]}, $row->[$column]);
		}
	}
	
	#print "retro $n";
	#print Dumper \@sorted_AoA;
	exit;
	
	return @sorted_AoA;
}


#######################
=head2 choose_transcript

 Usage   : $transcript = choose_transcript($transcript_choice, $transcripts);
 Function: Returns a single transcript from the list of transcripts of one gene ($transcripts), depending
           on the user's choice: either longest, shortest, or median length
 Returns : $transcript (GAL object)
 Args    : $transcript_choice (max|min|mid), $transcripts (GAL object, all transcripts of a gene)

=cut

sub choose_transcript {
	my ($transcript_choice, $transcripts) = @_;
	my $transcript;

	# User wants to analyze only shortest transcripts
	if ($transcript_choice eq 'min') {
		($transcript) = $transcripts->shortest;
	}
	
	# Default of User's choice median
	elsif ($transcript_choice eq 'mid') {
		my @transcripts = sort {$a->length <=> $b->length} $transcripts->all;
		my $middle = int(scalar(@transcripts)/2);
		$transcript = $transcripts[$middle];
	}
	
	# User want to analyze longest transcripts
	elsif ($transcript_choice eq 'max') {
		($transcript) = $transcripts->longest;
	}
	
	return $transcript;
}

#######################
=head2 print_batch_data

 Usage   : print_batch_data($FH, $type, $data_REF, @indices);
 Function: Prints mean or median (type) of each column of the given index list 
			from data_REF, into a tab separated line.
 Returns : -
 Args    : filehandle, 'mean' or 'median', reference to array of arrays (data columns),
			start index, end index

=cut

	sub print_batch_data {
		
		#Get filehandle, desired value type, data, list of columns
		my $FH			= shift @_;
		my $type		= shift @_;
		my $data_REF	= shift @_;
		my $indices_REF	= shift @_;
				
		my $last_i = ${$indices_REF}[-1];
		
		# Go through columns 

		foreach my $i (@{$indices_REF}) {
			my $value;
			
			if ($type eq 'mean') {
				# Get mean for each columns
				$value = my_mean(@{$data_REF}[$i]);
			}
			elsif ($type eq 'median') {
				# Get median for each columns
				$value = my_median(@{$data_REF}[$i]);
			}
		
			# And print it
			if ($value eq 'NA') {
				print $FH 'NA';
			} 
			else {
				printf $FH "%.5f", $value;
			}
			
			print $FH $t if !($i == $last_i);
		}
	}

#######################
=head2 calculate_CpG_oe

 Usage   : my $CpG_oe = calculate_CpG_oe($sequence);
 Function: Calculates the normalized CpG dinucleotide content [CpG observed/expected (o/e)] for
			the given sequence.
 Returns : $CpG_oe (floating integer)
 Args    : sequence of nucleotides

=cut
sub calculate_CpG_oe{
	# Written by Panagiotis Provataris
	
	my $seq = shift @_;
	
	# Count the number of occurrences for C, G and CG
		my $nr_C 		= $seq =~ s/C/C/gi;
		my $nr_G 		= $seq =~ s/G/G/gi;
		my $nr_CG		= $seq =~ s/CG/CG/gi;
		
		
	# Calculate frequencies of C, G and CG
		my $freq_C 		= $nr_C/length $seq;
		my $freq_G 		= $nr_G/length $seq;
		my $freq_CG		= $nr_CG/((length $seq) - 1) if ((length $seq) - 1) > 0;


	# Calculate CpG depletion
		my $CpG_dep 	= 'NA'; 
	
		if (($freq_C * $freq_G) > 0) {
			$CpG_dep	= sprintf "%.3f", ($freq_CG / ($freq_C * $freq_G));
		}
		
	return $CpG_dep;
}



########################################################################
#																	   #
#							FIN										   #
#																	   #
########################################################################

=head1 CONFIGURATION AND ENVIRONMENT


=head1 DEPENDENCIES

 COGNATE uses the following perl modules:
 5.010 strict warnings 
 Getopt::Long Data::Dumper
 FindBin FindBin::RealBin File::Path Cwd
 Bio::DB::Fasta 
 Number::Format List::Util Statistics::Basic
 GAL::Annotation
 GAL::List
  
 (Some) Modules in GAL/lib use the following further modules:

 Carp DBD::SQLite DBI Scalar::Util
 Set::IntSpan::Fast Statistics::Descriptive Text::RecordParser
 FileHandle IO::Prompt List::MoreUtils
 TAP::Harness Test::More Test::Pod::Coverage URI::Escape
 XML::LibXML::Reader


=head1 INCOMPATIBILITIES

 None reported.

=head1 BUGS AND LIMITATIONS

 No bugs have been reported, although there certainly are some.

 Please report any bugs or feature requests to <j.wilbrandt@zfmk.de>


=head1 AUTHOR

 Jeanne Wilbrandt  <j.wilbrandt@zfmk.de>

=head1 LICENSE AND COPYRIGHT

 This module is free software; you can redistribute it and/or
 modify it under the same terms as Perl itself (GPLv3).

=head1 DISCLAIMER OF WARRANTY

 BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
 FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
 OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
 PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
 EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
 ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
 YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
 NECESSARY SERVICING, REPAIR, OR CORRECTION.

 IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
 WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
 REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
 LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
 OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
 THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
 RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
 FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
 SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
 SUCH DAMAGES.
