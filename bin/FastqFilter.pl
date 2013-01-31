#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Verbose;
use Verbose::ProgressBar;

use Fastq::Parser '0.07';
use Fastq::Seq '0.09';


our $Version = '1.02';

=head1 NAME 

FastqFilter.pl

=head1 DESCRIPTION

Read Fastq file, gather some stats (length distributions,
 Nxx-values, base composition, quality distribution).
Modify sequences (trim lcs, mask low quality).
Report sequences based on some filters (length, ids, number of
 sequences).

=cut

=head1 CHANGELOG

=head 1.03 [Thomas Hackl 2013-01-28]

=over

=item [Feature] C<--phred-transform> to switch phred offsets between 33 and 64

=item [BugFix] Did set sliding window size instead of min value.

=item [Feature] IDS can have leading '@' or '>'.

=item [Change] C<--min-length/--max-length> filter are now executed
 after any trimming stuff.

=item [Feature] Autodetect <--phred-offset>.

=item [Feature] Sliding window based C<--trim>. 

=item [Renaming] C<--slice-lcs> to C<--trim-lcs>.

=back

=over

=item 1.02 [Thomas Hackl 2012-13-09]

Bugfix. Nx computation now correctly uses total length of filtered
 output, not total length of input.

=item 1.01 [Thomas Hackl 2012-11-09]

Bugfix. Renamed C<< reads_ids >> to C<< read_ids >>.

=item 1.00 [Thomas Hackl 2012-11-08]

Major Refactoring. 

=item 0.06 [Thomas Hackl 2012-11-06]

Bugfix. --fasta now accuratly also prints the descrition of the header line.

=item 0.05 [Thomas Hackl 2012-10-02]

Added --slice_lcs

=item 0.04 [Thomas Hackl 2012-09-27]

Removed --n-content, added --base-content, which allows you to specify 
 particular bases or group of base in a "," separated list to be scanned for.

=item 0.03 [Thomas Hackl 2012-08-31]

Added --ids, --ids-exclude

=item 0.02 [Thomas Hackl 2012-08-30]

Versioned, Verbose::ProgressBar, --min-length

=cut

=head1 TODO

=over

=item paired

=back

=cut



##------------------------------------------------------------------------##

=head1 SYNOPSIS

  perl FastqFilter.pl --in <FASTQFILE> <OPTIONS>
  
  cat <FASTQFILE> | perl FastqFilter.pl <OPTIONS>

=cut

=head1 OPTIONS

=cut

=over

=cut

my %opt;

=item [--in=<FASTA>]

Input Fastq file. Default STDIN.

=cut

$opt{'in=s'} = \(my $opt_in = undef);

=item [--out=<STRING>]

Output Fastq file. Default off. Specify '-' for STDOUT.

=cut

$opt{'out=s'} = \(my $opt_out = undef);

=item [--tsv=<STRING>/'-']

Print tab separated statistics of data AFTER filtering to pathname, 
 default off, specify '-' for STDOUT, but only if '--out' is not redirected 
 to STDOUT.

Column order is: number of sequences, total length, longest, shortest, 
 Nx-values, base-contents. No header is printed, so multiple runs can 
 simply be concatenated.

=cut 

$opt{'tsv=s'} = \(my $opt_tsv = '');

=item [--ids=<FILE>/--ids='-']

Pathname to file of sequence IDs or list of IDs to be reported. Reads
 comma-separated, whitespace and newline separated lists. Specify '-' to 
 read IDs form STDIN. Leading '>' or '@' are allowed but not necessary.

=cut 

$opt{'ids=s'} = \(my $opt_ids = '');

=item [--ids-exclude]

Exclude reads specified by --ids from being reported. Default off.

=cut 

$opt{'ids-exclude'} = \(my $opt_ids_exclude = undef);

=item [--min-length=<INT>]

Minimum sequence length.

=cut

$opt{'min-length=i'} = \(my $opt_min_length = 0);

=item [--max-length=<INT>]

Maximum length for sequences to be retrieved. Default off. Can be combined
 with max-count.
 
=cut 

$opt{'max-length=i'} = \(my $opt_max_length = undef);

=item [--fasta]

Output FASTA format. Default FASTQ.

=cut

$opt{'fasta'} = \(my $opt_fasta = undef);

=item [--phred-offset] [auto]

Phred offset for quality scores. Default auto-detect.

=cut

$opt{'phred-offset=i'} = \(my $opt_offset);

=item [--phred-transform]

Transform phreds from input offset to specified C<--phred-transform> offset,
 usually 33 to 64 or wise versa.

=cut

$opt{'phred-transform'} = \(my $opt_transform);

=item [--phred-mask]

Two values separated by ",", e.g "0,10" to mask all Nucleotides with
 phred below 10 with an "N".
 
NOTE: Requires C<--phred-offset>.

=cut

$opt{'phred-mask=s'} = \(my $opt_mask = undef);

=item [--trim=<INT1>,[<INT2>]]

Trim sequences to quality >= INT1 using a sliding window of size
 INT2, default 10. The sliding window allows to have a few positions
 below the <--trim> cutoff as long as the window mean is higher than
 INT1. It is made sure that a) below cutoff positions only occur
 within the remaining sequence, not at its start/end and b) windows
 never overlap eachother.

=cut

$opt{'trim=s'} = \(my $opt_trim = undef);

=item [--trim-lcs=<INT,INT,INT>]

Three values separated by ",", e.g. "30,40,50" to grep all stretches
 of quality >= 30 and minimum length 50 from the sequences.
Faster than C<--trim> yet breaks sequences even on a single low
 quality position.

NOTE: C<--trim> and C<--trim-lcs> cannot be combined.

=cut

$opt{'trim-lcs=s'} = \(my $opt_slice_lcs = undef);

=item [--N=<INT,INT...>]

Report Nx value (N50, N90...). Specify multiple values as comma separated
 STRING. Takes filter settings into account. Result might be ommited if to
 few sequences are retrieved. Default "50,90".

=cut 

$opt{'N=s'} = \(my $opt_N = '90,50');

=item [--base-content=<BASE(S),BASE(S),BASE(S),...>]

Count and print the relative amount of given bases. Takes a "," separated 
 list, each element of the list can be one or more bases. In the letter case
 the commulative amount of the individual bases is calculated.
 
Examples:
  --base-content=A,T,G,C,N   # content of regular bases
  --base-content=GC          # combined GC content

=cut

$opt{'base-content=s'} = \(my $opt_base_content = undef);

=item --[no]verbose

Toggle verbose, default on.

=cut 

$opt{'verbose!'} = \(my $opt_verbose = 1);

=item [--help]

Display this help

=cut

$opt{'help|?'} = \(my $opt_help);

=back

=cut


GetOptions(%opt) or pod2usage(1);
pod2usage(1) if $opt_help;
if($opt_max_length && $opt_min_length && ($opt_min_length > $opt_max_length)){
	pod2usage(exitval=>1, msg=>'min-length has to be smaller than max-length'); 
}elsif($opt_ids eq '-' && !$opt_in){
	pod2usage(exitval=>1, msg=>'Cannot read IDS and FASTA from STDIN');
}
if($opt_trim && $opt_slice_lcs){
	pod2usage(exitval=>1, msg=>'Cannot perform --trim and --trim-lcs at the same time'); 
}

##------------------------------------------------------------------------##

=head1 MAIN

=cut

my $V = Verbose->new(
	report_level => $opt_verbose || 0,
	line_width => 80
);


##------------------------------------------------------------------------##

=head2 check input file

=cut

my $fqp = Fastq::Parser->new(
	file => $opt_in, # defaults to STDIN if undef
)->check_format;
$V->exit($opt_in." does not look like FASTQ") unless $fqp;

$opt_offset = $fqp->guess_phred_offset unless $opt_offset;

pod2usage(msg=>'Could not guess phred offset, please specify', exitval => 1) unless $opt_offset;

$V->verbose('Detected FASTQ format, phred-offset '.$opt_offset);

my $VPB = Verbose::ProgressBar->new(
	size => $fqp->fh()
);


# output file
my $ofh;
if($opt_out){
	if($opt_out eq '-'){
		$ofh = \*STDOUT;
	}else{
		open($ofh, ">", $opt_out) || $V->exit("Can't open output read file: '$opt_out'");
	} 
}



##------------------------------------------------------------------------##

=head2 prepare som filter

=cut

my @opt_base_content = split(/\s*,\s*/, $opt_base_content) if $opt_base_content;
my @opt_N = sort{$a <=> $b}(split(/\s*,\s*/, $opt_N)) if $opt_N;
my %IDS = read_ids() if $opt_ids;

my $tsv;
if($opt_tsv){
	if($opt_tsv eq '-' && (!defined($opt_out) || $opt_out ne '-')){
		$tsv = \*STDOUT;
	}else{
		open($tsv, ">", $opt_tsv) || $V->exit("Can't open tsv stats file: '$opt_tsv'");
	} 
}

if($opt_mask){
	my($from, $to, $length) = split(',', $opt_mask);
	Fastq::Seq->Qual_low_range($from, $to, $opt_offset);
}

if($opt_slice_lcs){
	my($from, $to, $length) = split('\s*,\s*', $opt_slice_lcs);
	Fastq::Seq->Qual_lcs_range($from, $to, $opt_offset);
	Fastq::Seq->Qual_lcs_min_length($length);
}

if($opt_trim){
	my ($min, $wsize) = (split(',', $opt_trim), 10);
	Fastq::Seq->Qual_window_min($min);
	Fastq::Seq->Qual_window_size($wsize);
}

# TOTAL
my $total_count;	# total count
my $total_length = 0;
my $total_longest = 0;
my $total_shortest = undef;

# FILTERED
my @L; # lengths of seqs
my $filtered_length = 0;
my $filtered_longest = 0;
my $filtered_shortest = undef;
my %Ns; #Nx values
my %BCs; # base contents



##------------------------------------------------------------------------##

=head2 loop through file and apply filter

=cut


# Fasta parser
$V->verbose("Reading Input from ".( $opt_in ? $opt_in : "STDIN" ));

# loop through FASTQ
while (my $fq = $fqp->next_seq()){
	$total_count++;
	$VPB->update();

	my $length = length($fq->seq);
	
	# TOTAL STATS
	# length stuff
	$total_length+=$length;
	$total_longest = $length if $length > $total_longest;
	$total_shortest = $length if (!defined ($total_shortest) || $length < $total_shortest);
	
	
	# FILTER
	# ids
	if(keys %IDS){
		if ($opt_ids_exclude){
			next if $IDS{$fq->id}
		}else{
			next unless $IDS{$fq->id}
		}
	}
	
	# slice by lcs
	my @fq;
	if($opt_trim){
		@fq = $fq->slice_seq($fq->qual_window(1));
	}elsif($opt_slice_lcs){
		@fq = $fq->slice_qual_lcs(1);
	}else{
		@fq = ($fq);
	}
	
	foreach my $fq (@fq){
		my $length = length($fq->seq);
		# min-length/max-length
		if($opt_min_length){ next if $length < $opt_min_length }; 
		if($opt_max_length){ next if $length > $opt_max_length };
		
		# mask low quality
		$opt_mask && $fq->mask_qual_low();

		# OUTPUT
		if($opt_out){
			if($opt_fasta){
				print $ofh sprintf(">%s\n%s\n", substr ($fq->seq_head, 1), $fq->seq)
			}elsif($opt_transform){
				print $ofh $fq->phred_transform()->string();
			}else{
				print $ofh $fq->string();
			}
		}
		
		# OUTPUT STATS
		# store length for Nx computation
		push @L, $length;
		# length stuff
		$filtered_length+=$length;
		$filtered_longest = $length if $length > $filtered_longest;
		$filtered_shortest = $length if (!defined ($filtered_shortest) || $length < $filtered_shortest);
		# base content
		if(@opt_base_content){
			foreach (@opt_base_content){
				$BCs{$_} += $fq->base_content($_);
			}
		}
	}
}

$VPB->finish();


=head2 print summary

=cut

$V->hline();
$V->verbose('Input');
$V->verbose(sprintf("%-15s %10d #","Sequences", $total_count));
$V->verbose(sprintf("%-15s %10d bp", "Total", $total_length));
$V->verbose(sprintf("%-15s %10d bp", "Longest", $total_longest));
$V->verbose(sprintf("%-15s %10d bp", "Shortest",$total_shortest));
$V->hline();


unless(keys @L){
	$V->verbose("No sequences matching filters found");
	$V->exit("\n");
}

#Nx computation
if(@opt_N){
	@L = sort{$b <=>$a}@L;
	my @ls = map{$filtered_length * ($_/100)}@opt_N;
	
	my $lc = 0;
	foreach my $l (@L){
		$lc+= $l;
		if($lc >= $ls[0]){
			shift @ls;
			$Ns{shift @opt_N} = $l;
			last unless @ls;
		}
	}
}


$V->verbose("Filtered");
$V->verbose(sprintf("%-15s %10d #","Sequences", scalar @L));
$V->verbose(sprintf("%-15s %10d bp", "Total", $filtered_length));
$V->verbose(sprintf("%-15s %10d bp", "Longest", $filtered_longest));
$V->verbose(sprintf("%-15s %10d bp", "Shortest", $filtered_shortest));

my @Ns;
foreach (sort{$a <=> $b}keys %Ns){
	$V->verbose(sprintf("%-15s %10d bp", "N$_", $Ns{$_}));
	push @Ns, $Ns{$_};
}

# base content
my @BCs;
if(keys %BCs){
	foreach my $k(sort keys %BCs){
		$V->verbose(sprintf("%-15s %10d bp %10.2f %%", "[$k]", $BCs{$k}, (100*$BCs{$k}/$total_length)));
		push @BCs, $BCs{$k};
	}
}


$V->hline();
$V->nline();



if($tsv){
	print $tsv join("\t", scalar @L, $filtered_length, $filtered_longest, $filtered_shortest, @Ns, @BCs),"\n";
}



##------------------------------------------------------------------------##

=head1 Methods

=cut

=head2 read_ids

=cut

sub read_ids{
	my @IDS;
	if($opt_ids eq '-'){
		while(<STDIN>){
			push @IDS, split(/[\s,\n]+/, $_);
		}
	}else{
		open(IDS, $opt_ids) || $V->exit("Can't open output ids file: '$opt_ids'");
		while(<IDS>){
			push @IDS, split(/[\s,\n]+/, $_);
		}
		close IDS;
	}
	map{s/^[>@]//}@IDS;
	my %IDS;
	@IDS{@IDS}=map{1}@IDS;
	return %IDS;
}



=head2 chars2phred(chars => '', [offset=>64])

Convert ascii quality sequence characters to integer phred values according 
 to specified phred offset. Returns STRING if chars where given as STRING or
 ARRAY if chars where given as ARRAY. The input order is kept in both cases.

PARAMETER

=over

=item chars

STRING or ARRAYREF

=item offset [64]

INT

=back

=cut

sub chars2phred{
	my %opt = (
		offset => 64,
		chars => '',
		@_
	);
	
	my $string_flag;
	
	unless (ref $opt{chars} eq 'ARRAY'){
		$string_flag++;
		$opt{chars} = [split(//, $opt{chars})];
	}
	
	my @phreds = map{ ord($_)-$opt{offset} }@{$opt{chars}};
	
	return $string_flag ? join('', @phreds) : @phreds;
};


=head2 phred2base_call_accuracy(phred)

Returns the base call accuracy probability for a given phred score.

=cut

sub phred2base_call_accuracy{
	my $phred = shift;
	return 1-(10**(-$phred /10));
}

=head2 phred2char(phred score, [offset])

Returns ascii char for given phred score based on offset, default 64;

=cut

sub phred2char{
	my ($phred, $offset) = (@_, 64);
	return chr($phred+$offset);
}

