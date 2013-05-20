#!/usr/bin/env perl

# $Id: ChimeraToSeqFilter 644 2013-01-18 16:04:44Z dumps $

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

use Verbose;
use Verbose::ProgressBar;

use Data::Dumper;
$Data::Dumper::SortKeys=1;
$Data::Dumper::SortKeys=1;

our $VERSION = '0.01';

=head1 NAME

ChimeraToSeqFilter - convert chimera annotations to SeqFilter --substr 
 parameter

=cut

=head1 CHANGELOG

=head1 0.01

=over

=item [Init]

=back

=cut

=head1 TODO

=over

=back

=cut



##------------------------------------------------------------------------##

=head1 SYNOPSIS

  perl ChimeraToSeqFilter --in <*.chimera.tsv> <OPTIONS> --out <OUTFILE>
  
=cut

=head1 OPTIONS

=cut

=over

=cut

my %opt;

=item --in=<TSV>

Input proovread chimera tsv file. Required.

=cut

$opt{'in=s'} = \(my $opt_in = undef);

=item [--out=<STRING>]

Output file. Default STDOUT.

=cut

$opt{'out=s'} = \(my $opt_out = undef);

=item [--min-score=<INT>=0.01]

Minimum score of a chimera breakpoint to be accepted.

=cut

$opt{'min-score=s'} = \(my $opt_min_score = 0.01);

=item [--deviation=<INT>=20]

The number of breakpoint flanking that should be trimmed.

=cut 

$opt{'trim-length=i'} = \(my $opt_trim = 20);

=item --verbose=<INT>

Toggle verbose level, default 2, which outputs statistics and progress.
Set 1 for statistics only or 0 for no verbose output.

=cut 

$opt{'verbose=i'} = \(my $opt_verbose = 2);

=item --quiet

Omit all verbose messages. The same as --verbose=0, Superceeds --verbose settings.

=cut

$opt{'quiet'} = \(my $opt_quiet);

=item [--help]

Display this help

=cut

$opt{'help|?'} = \(my $opt_help);

=back

=cut



GetOptions(%opt) or pod2usage(1);

# First arg is in
if(@ARGV && !($ARGV[0] =~ /^-/) && !$opt_in){
	$opt_in = $ARGV[0];
}

my $oh = \*STDOUT;
$opt_out && open($oh, '>', $opt_out) or die "$!: $opt_out";
select $oh;

pod2usage(1) if $opt_help;
$opt_in || pod2usage(exitval=>1, msg=>'Input file required');
$opt_verbose = 0 if $opt_quiet; 



##------------------------------------------------------------------------##

=head1 MAIN

=cut

my $V = Verbose->new(
	level => 1,
	report_level => $opt_verbose,
	line_width => 80,
	
);

##------------------------------------------------------------------------##

=head2 check input file

=cut

$V->verbose("Reading Input from ".( $opt_in ? $opt_in : "STDIN" ));

my %read = (id => '', coords =>[]);

my $bpc;

open(CHIM, '<', $opt_in) or $V->exit($!);
my $header = scalar <CHIM>;


while(<CHIM>){
	chomp();
	my ($id, $from, $to, $score) = split(/\t/, $_);
	if($id ne $read{id}){# new read
		if(@{$read{coords}}){
			my @coords = @{$read{coords}};
			unshift(@coords, 0);
			my $i;
			for($i=0; $i<@coords-1; $i+=2){
				printf ("%s\t%s\t%s\n", $read{id}, $coords[$i], $coords[$i+1]);
			}
			printf ("%s\t%s\n", $read{id}, $coords[$i]);
		}
		# create new read
		%read = (id => $id, coords => []);

	}else{
		if($score >= $opt_min_score){
			push @{$read{coords}}, $from, $to;
		}
	}

}

close CHIM;





=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut

