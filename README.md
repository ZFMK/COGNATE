******************** COGNATE README *********************

A. Introduction

B. COGNATE setup

C. COGNATE options

D. COGNATE input

E. COGNATE output (file indices)

F. COGNATE setup test with sample data

G. Output usage examples (in progress)

H. Notes on GFF3 files

I. References


COGNATE v1.0	February 2017

Contact:	j.wilbrandt@leibniz-zfmk.de

***********************************************************

A. Introduction
---------------

COGNATE (COmparative GeNe AnnoTation charactErizer) is a tool to simultaneously analyze 
a given protein-coding gene annotation and the corresponding assembled sequences of a genome.
It allows a quick and easy extraction of basic genome feature and gene repertoire data, 
such as exon and intron lengths, counts, as well as the coding amount (i.e., the added length 
of all coding sequences (CDSs)). Thus, COGNATE is a tool to primarily describe a genome and its 
annotation of protein-coding genes, which is an essential prerequisite for comparative and meta-analyses
of genome and gene structure. 

COGNATE infers the following main parameters:
 - summary counts of the analyzed features
 - counts of single CDS/exon genes, strand-mix genes (CDS/exon/intron on different strand than 
   transcript), alternative spliceforms per transcript
 - length statistics (nucleotide sequences/amino acid sequences), including N50/L50, 75/L75, N90/L90
 - intron length distribution as suggested by [1]
 - percental GC content statistics in two different ways, namely 
        - using a calculation that explicitly considers IUPAC ambiguity codes 
	        (G, C, S per total length excluding N, R, Y, K, M, B, D, H, V)
	      - using the previously prevailing calculation of GC per total length, which is inappropriate 
          for genome comparisons due to its dependence on assembly quality
 - statistics of CpG dinucleotide depletion (CpG observed/expected), normalized by C and G content 
   of the respective region [2]
 - density statistics (ratio of the length of a feature covered by another, number-wise)
 - coverage statistics (ratio of the length of a feature covered by another, length-wise).

The output design is focused on clarity and the combination of overview display and detailed 
parameter evaluation to make COGNATE as useful as possible within its designated scope.

Usage of the tool COGNATE requires a working perl installation and
GAL (Genome Annotation Library, https://github.com/The-Sequence-Ontology/GAL), 
which is originally provided by the Sequence Ontology [3], written by Barry Moore. 
The GAL/lib/ directory is included with minor changes as GAL/ in the 
distribution of COGNATE and thus requires no further installation.


***********************************************************

B. COGNATE setup
----------------

The COGNATE distribution is released as a compressed archive file (COGNATE.tar.gz) for
download. 

Extracting the files to your current directory `tar -zxvf COGNATE.tar.gz` will
create the directory COGNATE/, containing the required files.

The used GAL libraries provided with the package need no further installation:
 - GAL::Annotation
 - GAL::List
  
Before you begin, you will need to make sure that the following required perl modules/libraries
are installed:

 - FindBin
 - Bio::DB::Fasta
 - Number::Format
 - List::Util
 - Statistics::Basic
 - Statistics::Descriptive
 - DBD::SQLite
 - DBI
 - DBIx::Class
 - Scalar::Util
 - Set::IntSpan::Fast
 - Statistics::Descriptive
 - Text::RecordParser
 - IO::Prompt
 - List::MoreUtils
 - TAP::Harness
 - Test::More
 - Test::Pod::Coverage
 - URI::Escape
 - XML::LibXML::Reader

You can do this by running COGNATE_CheckDep.pl: `$ perl COGNATE_CheckDep.pl`

If there are modules missing, you need to install them, for example via cpanminus. If this is
not installed yet, type in your terminal: `$ sudo apt-get install cpanminus`

If you have cpanminus installed, type for each MODULE to install it: `$ sudo cpanm MODULE`


***********************************************************

C. COGNATE arguments and options
--------------------------------

### 1 - Mandatory arguments 
```
 --gff /PATH/TO/FILE.gff3	use gff3 file with absolute path as input
 --fasta /PATH/TO/FILE.fa	use fasta file with absolute path (accepted suffixes: .fa,.fn,.fasta,.fas,.fna)
 ```
 
 alternatively:
 ```
 --input /PATH/TO/COGNATE.input	use COGNATE.input file (specifics see below)
 ```
 
### 2 - Optional arguments 

```
-n|--name STRING		Name of the genome/species, used to name the output (e.g, my_species). 
				 Default: genome_# (number increments automatically)
-b|--batch STRING		Write summary line to batch output file, STRING is the name of the analyzed
				 set (e.g., my_genome_set).
-w|--workingdir	/PATH/		Determine working directory.
				 Default: current directory.
-o|--overwrite			Overwrite existing COGNATE directories (for NAME), no questions asked.
-f|--no_filecheck		Switch off validity check of input files.
-l|--transcript_length STRING	Evaluate only the max=longest, min=shortest or mid=middle transcript.
				 Default: max.
-t|--count_table		Switch on the overview table (printed in STDOUT).
--print LIST			Comma-separated list (e.g., 1,2,3) of output file indices (see below), which will be printed
				 Default: all files are printed.
--dont_print LIST		Comma-separated list (e.g., 1,2,3) of output file indices (see below), which will not be printed
-h|--help 			Print help (more info: perldoc COGNATE_v1.0.pl).
```

***********************************************************

D. COGNATE input
------------------

Files:
- GFF (feature file/annotation) in gff3 format.
- FASTA (genome sequences) in fasta format. Make sure that the sequences are not in one line but
   include line-breaks (after, e.g., 80 bases). All sequence lines must have the same length except the
   last one.

File supply:
a) If the option `--gff` and `--fasta` (additional recommendation: --name) are chosen, both files need to
   be provided with absolute paths.
 
b) As an alternative and to facilitate batch submissions of multiple genomes, you can give COGNATE a
   input file using the option `--input`. This file needs to follow these rules:
   - file name end with COGNATE.input
   - lines beginning with # are ignored (comments)
   - each line contains input data for one genome/species in three tab-separated columns (see as an
     example the file triple-dmel-ex_COGNATE.input in the Example_data dir):
     * name
     * /PATH/TO/FILE.gff3
     * /PATH/TO/FILE.fa
     * Example content:
     ```
     # NAME		GFF3			FASTA
     species_1	/PATH/annotation_1.gff3	/PATH/genome_1.fa
     genome_b	/PATH/annotation_b.gff3	/home/user/Desktop/genome_b.fa
     ```
     
***********************************************************


E. COGNATE output
------------------

COGNATE produces per default 21 output files. Their production can be controlled using
either the option `--print` or `--dont_print`. The argument for this option has to be a
comma-separated, white-space-free list of numbers (= indices of output files). Find 
the file indices in the lists below.

Successful execution of COGNATE will create a directory named COGNATE_NAME, where NAME is your 
assigned name for the genome-specific run. BATCH is the string given to `--batch`. Only with this
option, files 14-20 will be created!
The newly created output directory (1) and the working directory (2) will contain these files:

### 1 - Files in COGNATE_NAME 
`# Index		Name & content`

 00	NAME_00_analyzed_transcripts.fa
 
 01	NAME_01_summary.tsv
 
 02	NAME_02_scaffold_general.tsv
 
 03	NAME_03_scaffold_transcripts.tsv
 
 04	NAME_04_scaffold_CDSs.tsv
 
 05	NAME_05_scaffold_exons.tsv
 
 06	NAME_06_scaffold_introns.tsv
 
 07	NAME_07_transcript_general.tsv
 
 08	NAME_08_transcript_CDSs.tsv
 
 09	NAME_09_transcript_exons.tsv
 
 10	NAME_10_transcript_introns.tsv
 
 11	NAME_11_CDSs.tsv
 
 12	NAME_12_exons.tsv
 
 13	NAME_13_introns.tsv
 


### 2 - Files in the working directory (not affected by --overwrite) 
`# Index		Name`

 14	BATCH_14_batch_general.tsv
 
 15	BATCH_15_batch_scaffold-means.tsv
 
 16	BATCH_16_batch_scaffold-medians.tsv
 
 17	BATCH_17_batch_transcript-means.tsv
 
 18	BATCH_18_batch_transcript-medians.tsv
 
 19	BATCH_19_batch_component_sizes.tsv
 
 20	BATCH_20_batch_bash-commands.txt
 

More than 280 parameters will be evaluated and summarized. Which parameter is 
contained in which file can be found in Table S1 of the publication.


### Notes 

 - Of most parameters both mean and median are given. The mean is often inappropriate due to non-normal distribution
   of data. Thus, for comparisons, the median is recommended.
 - A median of means or medians for a certain parameters results from the calculation of the median/mean values per 
   structure entity, which in turn were calculated for all sub-structures. For example: CDS length for one CDS -> 
   median CDS length per transcript for one transcript -> median of median CDS length per transcript for the whole annotation.

 - For GC content, two types are calculated: total (gc/length) and (non-)ambiguity (gcs/length-NRYKMBDHV) [noAm].
   The latter GC content is not dependent on assembly quality. Both types include softmasked sequences.
   Thus, 'GC content Without ambiguity' means that the length of a nucleotide sequence was tallied excluding the bases 
   N, R, Y, K, M, B, D, H, V (IUPAC codes for all bases that are not G or C (S = G/C)) and the GC content
   calculated as count of G, C, S / length without ambiguity.

 - Protein length is calculated without stop codon (*).

 - COGNATE only evaluates one transcript per gene, thus no isoforms although they are given in NCBI gffs. Alternative 
   splice-forms are counted, though.
 - Count of alternative spliceforms (files 01, 07) = count of annotated transcripts per gene minus 1.

 - L90[genes] (file 01) is the count of largest scaffolds required to find at least 90 % of the annotated and by COGNATE 
   analyzed genes.

 - CpG o/e is the ratio of CpG-dinucleotide depletion and calculated as (frequency of CG/(frequence of C * frequency of G)), 
   where the frequency is the count of a (di)nucleotide in a sequence / length of this sequence.

 - Count of strand-mix genes (file 01) = count of transcripts where a CDS/exon/intron differs from the transcripts strand. 
 - Strands of transcripts (file 07) are given as +/-. Individual strandedness for CDSs/exons/introns can be found 
   in the respective file (11/12/13) and is given as +/-, in case of a conflict with the transcript strandedness as +!/-!.
 - Transcripts on + : - strand (file 03) = percentage of all transcripts on + strand and percentage of all transcripts 
   on - strand (normalization against total count of non-isoform transcripts).

 - COGNATE checks for overlapping transcripts (same SCS, same strand, coordinates overlap). If there are any, the respective
   section is appended to the summary (file 01).


***********************************************************

F. COGNATE setup test with sample data
--------------------------------------

Sample data are provided to test your COGNATE setup. Execute the following steps 
and compare the final output with the provided files in /Example_output/.

 1. Change directory to unpacked COGNATE package

 	`$ cd COGNATE/`
  
 2. Adjust /PATH/s according to your situation (also inside triple-dmel-ex_COGNATE.input for run option b).

 3. Run COGNATE on the sequence and annotation files provided in /Example_data/. 
    Use either two explicit input files for a single-genome run (a) or the COGNATE.input file for a
    multi-genome run (b).
	
	a) `$ perl COGNATE_v1.0.pl --gff /PATH/Example_data/Dmel_example_annotation.gff --fasta /PATH/Example_data/Dmel_example_scaffolds.fa --name dmel_ex --batch Example`
  
 	b) `$ perl COGNATE_v1.0.pl --input /PATH/Example_data/triple-dmel-ex_COGNATE.input --batch Example`
  
 3. Compare the output with the provided example output in /Example_results/.


### Notes 

 - If you want to view the .tsv result files using LibreOffice (or similar tools), be sure to
   * use Unicode (UTF-8) character set
   * check 'Tab' at Separator Options
   * do not check 'Merge delimiters'
   * do not check 'Comma separated'

 - To interpret the output terms, reading the publication and the glossary (Table S2) included therein 
   may be hepful.


***********************************************************

G. Output usage examples (in progress)
--------------------------------------

 - concatenate columns of interest
   * concatenate all columns of one kind for all species of the batch
   * sort columns of interest
   * --> histograms

 - mini R tutorial
	 * read in data
	 * plot data
	 * name axes

***********************************************************

H. Notes on GFF3 files
----------------------
Although COGNATE offers to check the input gff3 for validity, it is always a good idea to know your input.

The first line of a proper GFF3 file has to look like this:
`##gff-version 3`

For more information on GFF(3) properties, have a look here: 
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

Some annotation pipelines may include organellar genes into the annotation or declare pseudogenes as such only 
in column 9 (e.g., 'pseudo=true'). Here, CanonGff3 (http://aegean.readthedocs.io/en/v0.16.0/canon.html) can be helpful.


***********************************************************

I. References
-------------
1. Elango N, Hunt BG, Goodisman MAD, Yi SV. DNA methylation is widespread and associated with 
	differential gene expression in castes of the honeybee, Apis mellifera. PNAS. 2009;106:11206–11. 
2. Roy SW, Penny D. Intron length distributions and gene prediction. Nucl. Acids Res. 2007;35:4737–42. 
3. Eilbeck K, Lewis SE, Mungall CJ, Yandell M, Stein L, Durbin R, et al. The Sequence Ontology: a tool 
	for the unification of genome annotations. Genome Biology. 2005;6:R44. 





***********************************************************
www.zfmk.de/en/cognate


