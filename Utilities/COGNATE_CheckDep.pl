#!/usr/bin/perl


use strict;
use warnings;

my $n = "\n";

########### Check COGNATE dependencies #################################
# - require modules needed by COGNATE
# - print either OK if present
# - or MISSING: if not.
# -- collect MISSING and print together at the end
# -- suggest cpanm
########################################################################

my @deps = qw(FindBin 
				Bio::DB::Fasta
				Number::Format
				List::Util 
				Statistics::Basic
				Statistics::Descriptive
				DBD::SQLite
				DBI
				DBIx::Class
				Scalar::Util
				Set::IntSpan::Fast
				Statistics::Descriptive
				Text::RecordParser
				Text::Table
				IO::Prompt
				List::MoreUtils
				TAP::Harness
				Test::More
				Test::Pod::Coverage
				URI::Escape
				XML::LibXML::Reader
				Getopt::Long
				Data::Dumper
				FindBin
				FindBin::RealBin
				File::Path
				Cwd
				Carp
				FileHandle
				);

print "# Checking COGNATE dependencies...$n";

my @missing_deps;

foreach my $dep (sort @deps) {
	if  (eval {"require $dep";1;}) {
		print "OK: $dep$n";
	}
	else {
		push @missing_deps, $dep;
	}
}

if (scalar @missing_deps > 0) {
	print "MISSING modules:$n";
	print join ($n, sort @missing_deps), $n;
	print "${n}Install these using for example 'sudo cpanm MODULE'. See COGNATE readme for further details.";
}


print "# Done$n";
