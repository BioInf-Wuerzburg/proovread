package Sam::Seq;

use warnings;
use strict;

# $Id$

use List::Util;

# preference libs in same folder over @INC
use lib '../';

use Verbose;

use Sam::Parser;
use Sam::Alignment ':flags';

use Fastq::Seq;
use Fasta::Seq;


our $VERSION = '0.09';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;


=head1 NAME 

Sam::Seq.pm

=head1 DESCRIPTION

Class for handling sam reference sequences and its aligned reads.

=cut

=head1 SYNOPSIS

  use Sam::Seq;
  
=cut

=head1 CHANGELOG

=head2 0.09

=over

=item [Feature] Chimera detection

=item [Change] Replaced C<each> and C<keys> statements on reference
 values, which are considered experimental features and would 
 require Perl 5.014. Also removed C<use 5.014>. 

=back

=head2 0.08

=over

=item [Change] In order to implement InDelTaboo correctly, the states
 of an alignment are now computed independently of the state matrix, then 
 trimmed and after that added to the global matrix.

=item [Feature] C<< $Sam::Seq::InDelTaboo >> and C<< Sam::Seq->InDelTaboo() >>

=item [Feature] length normalize score to handle varying short read lengths

=back

=head2 0.07

=over

=item [Change] Added svn:keywords

=back

=head2 0.06

=over

=item [Change]  C<< use lib '../' >>: Prioritise libs in same path over default
 @INC libs

=item [Feature] Added C<< Sam::Seq->BinSize >>, C<< Sam::Seq->BinMaxCoverage >>
 and C<< Sam::Seq->PhredOffset >> to controll the corresponding class attributes.

=back

=over

=item 0.05 [Thomas Hackl 2012-11-29]

New phred score calculation, using square root of coverage 
 probability of consensus state C<< freq2phred >>.

Ignore long insertions (length > 3) in query, since they are more
 likely to be mapping artefacts caused by cheap gap costs, than 
 actual pacbio sequencing errors (would be 4 or more insertions in
 a row...)

Merging of previously corrected FASTQs to _state_matrix now 
 uses C<< phred2freq >> to recompute frequency of chosen state.

=item 0.04 [Thomas Hackl 2012-11-02]

Set a minimum phred score of 10 (10% error rate) in case anything has been 
 aligned. 

=item 0.03

Reference and consensus sequences are now kept as Fasta::Seq/Fastq:Seq 
 objects within the Sam::Seq object;

=item 0.02

Merged Sam::Consensus functionality into this module.

Alignments are now stored in globally, local bins for local coverage 
 estimation only hold pointer to the global store.

Added Object ATTRIBUTES C<bin_size and bin_max_coverage>, which default to
 Class ATTRIBUTES C<$BinSize and $BinMaxCoverage>, respectively.

=item 0.01

Initial Alignment module. Provides Constructor, generic accessor 
 methods.

=back

=cut


=head1 TODO


=over

=item BUG

Use of uninitialized value $idx in hash element at /storage/genomics/scripts/lib\
/Sam/Seq.pm line 830, <$fh> line 54918.

=item Tests

=item Coverage

Currently there is a stand-alone coverage() methods, which uses the 
 _state_matrix and returns an array of absolute coverage values and the
 consensus methods stores a coverage in phred style (chr+offset) in 
 $self->cov seq. This is redundant...



=back

=cut

##------------------------------------------------------------------------##


=head1 Class ATTRIBUTES

=cut

=head2 $V

Verbose messages are handled using the Verbose.pm module. To 
 customize verbose message behaviour, overwrite the attribute with
 another Verbose object created with the Verbose module.

=cut

our $V = Verbose->new();

=head2 $BinSize [20]

=cut

our $BinSize = 20;

=head2 $MaxCoverage [50]

=cut

our $MaxCoverage = 50;

=head2 $PhredOffset [33]

=cut

our $PhredOffset = 33;

=head2 InDelTaboo [0.1]

Trim reads to prevent insertions/deletions within the first/last 
 InDelTaboo fraction of the read. N=0 deactivates the feature.

=cut

our $InDelTaboo = 0.1;

=head2 $Trim

Boolean. Deactivate trimming completely, including leading/trailing indels 
 and InDelTaboo.

=cut

our $Trim = 1;


=head2 %Freqs2Phreds

=cut

our %Freqs2Phreds;
@Freqs2Phreds{0..31} = (map{int((($_/50)**(1/2)*50)+0.5)}(0..39));
@Freqs2Phreds{32..100} = (40)x69;


our %Phreds2Freqs = reverse %Freqs2Phreds;

=head1 Class METHODS

=cut

=head2 BinSize

Get/Set $Sam::Seq::BinSize. Default 20.

=cut

sub BinSize{
	my ($class, $size) = @_;
	$BinSize = $size if defined $size;
	return $BinSize;
}

=head2 MaxCoverage

Get/Set $Sam::Seq::MaxCoverage. Default 50.

=cut

sub MaxCoverage{
	my ($class, $cov) = @_;
	$MaxCoverage = $cov if defined $cov;
	return $MaxCoverage;
}

=head2 PhredOffset

Get/Set $Sam::Seq::PhredOffset. Default 33.

=cut

sub PhredOffset{
	my ($class, $offset) = @_;
	$PhredOffset = $offset if defined $offset;
	return $PhredOffset;
}

=head2 InDelTaboo

Get/Set $Sam::Seq::InDelTaboo. Default 10.

=cut

sub InDelTaboo{
	my ($class, $indeltaboo) = @_;
	$InDelTaboo = $indeltaboo if defined $indeltaboo;
	return $InDelTaboo;
}

=head2 Trim

Get/Set $Sam::Seq::Trim. Default TRUE.

=cut

sub Trim{
	my ($class, $trim) = @_;
	$Trim = $trim if defined $trim;
	return $Trim;
}


=head2 Freqs2Phreds

Convert a LIST of frequencies to a LIST of phreds based on precomputed values
 from 

=cut

sub Freqs2Phreds{
	my $class = shift;
	return map{
		$_ = 100 if $_> 100;
		$Freqs2Phreds{int($_)};
	}@_;
}


=head2 Phreds2Freqs

Convert a LIST of Phreds to a LIST of frequencies based on precomputed values.

=cut

sub Phreds2Freqs{
	my $class = shift;
	return map{
		$_ = 40 if $_> 40;
		$Phreds2Freqs{int($_)};
	}@_;
}


##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

Create a Sam::Seq object. Takes key => value representation. 
Returns a Sam::Seq object.

  ## defaults
  id => undef,		
  len => undef,           # length of the reference sequence, required			
  ref => undef,           # reference seq object, Fasta::Seq or Fastq::Seq
  con => undef,           # consensus seq, Fastq::Seq (+cov)
  ref_merge_regions => [],# regions in the ref seq, to be included in 
                          #  consensus calling, ARRAY of Tuples (start, offset)
  max_coverage => 50,     # assumed maximum coverage value for consensus 
                          #  quality value calculation
  sam => undef,           # associate with a Sam::Parser object, only an
                          #  index positions to this file will be stored  
                          #  in _aln, data is assumed to be stored in this
                          #  file at respective position and will be 
                          #  retieved by reading from the Sam::Parser at 
                          #  this position
                          # default no file -> everything in ram, _aln 
                          #  contains real Sam::Alignment objects.
  is => undef,            # see Sam::Parser->is()
                          
  bin_size => $Sam::Seq::BinSize,
  bin_max_coverage => $Sam::Seq::BinMaxCoverage,
  phred_offset => $Sam::Seq::PhredOffset,
  

=cut

sub new{
	my $class = shift;
	my $self;
	
	$self = {
		## defaults
		id => undef,		
		len => undef,			# length of the reference sequence, required			
		ref => undef,			# reference seq object, Fasta::Seq or Fastq::Seq
		con => undef,			# consensus seq, Fastq::Seq (+cov)
		ref_con_regions => [],  # regions in the ref seq, to be included in 
								#  consensus calling, ARRAY of Tuples (start, offset)
		sam => undef,			# associate with a Sam::Parser object, only an
								#  index positions to this file will be stored  
								#  in _aln, data is assumed to be stored in this
								#  file at respective position and will be 
								#  retieved by reading from the Sam::Parser at 
								#  this position
								# default no file -> everything in ram, _aln 
								#  contains real Sam::Alignment objects.
		max_coverage => $MaxCoverage,
		bin_size => $BinSize,
		bin_max_bases => $BinSize * $MaxCoverage,
		phred_offset => $PhredOffset,
		max_coverage => 50,
		is => undef,
		## custom overwrites
		@_,
		## protected				
		_is => sub{return 1},	# don't set directly, use is()
		_alns => {},			# reads or idx pos of reads in sam
		_bin_alns => undef,		# position bins containing aln refs
		_bin_scores => undef,	# position bins containing nscores of alns in _bin_alns
		_bin_bases => undef,
		_bin_lengths => undef,
		_aln_idc => 0,			#
		_state_matrix => [],	# consensus state matrix
		_states => {			# consensus states
			A => 0,
			T => 1,
			G => 2,
			C => 3,
			'-' => 4,
			N => 5,
			# .. complex states, dynamically added
		}
		
	};


	bless $self, $class;

	# prepare is test routine
	$self->is($self->{is}) if $self->{is};#
	
	# init bins
	$self->_init_read_bins;

	return $self;
}







##------------------------------------------------------------------------##

=head1 Public METHODS

=cut

=head2 add_aln_by_score

Add a sam alignment object of a mapped illumina read to the object, based 
 on position and score. Takes a file position as second object. If provided
 this position is stored instead of the actual alignment, the position is 
 assumed to be a index pointing to the alignment in the file indicated by
 C<< $sam_seq->sam >>. Returns internal alignment id (>0) if aln has been 
 added, else 0.

  $sam_seq->add_aln_by_score($aln);
  $sam_seq->add_aln_by_score($aln, tell($sam_fh);


=cut

sub add_aln_by_score{
	my ($self, $aln) = @_;
	
	my $bin = $self->bin($aln);
	my $bases = length($aln->seq);
	my $nscore = $aln->opt('AS') / $bases;
	
	# if bin_bases are full, check if new nscore is good enough
	if( $self->{_bin_bases}[$bin] > $self->{bin_max_bases} ){
		# ignore scores, that are too low
		if( $nscore <= $self->{_bin_scores}[$bin][-1] ){
			return 0; 
		}else{ # sufficient score
			# remove lowest scoring aln before adding new one
			#  from score bins
			pop(@{$self->{_bin_scores}[$bin]});
			#  from length bin
			my $rm_bases = pop(@{$self->{_bin_lengths}[$bin]});
			
			# remove aln from global store 
			delete($self->{_alns}{
					# and remove aln id from aln bins
					pop(@{$self->{_bin_alns}[$bin]})
				});
				
			# adjust bin_bases by length length difference of new and old alignment
			$self->{_bin_bases}[$bin] += ($bases - $rm_bases);
		}
	}else{
		# new aln added w/o removal other aln -> add length to bin bases
		$self->{_bin_bases}[$bin] += $bases;
	}

	my $id = $self->add_aln($aln);
	
	# set score/id at the right place
	my $i = @{$self->{_bin_scores}[$bin]} - 1;
	$i-- while $i >= 0 && $nscore > $self->{_bin_scores}[$bin][$i];
	# store new  score and _id of aln at correct position
	splice(@{$self->{_bin_scores}[$bin]}, $i+1, 0, $nscore);
	splice(@{$self->{_bin_alns}[$bin]}, $i+1, 0, $id);
	splice(@{$self->{_bin_lengths}[$bin]}, $i+1, 0, $bases);

	return $id;
}

=head2 add_aln

Add a Sam::Alignment to Sam::Seq. If Sam::Seq->sam is set, it writes
 and stores the file position, else it stores the actual Sam::Alignment.

=cut


sub add_aln{
	my ($self, $aln) = @_;
	
	$self->{_alns}{++$self->{_aln_idc}} = $self->{sam} 
		? $self->{sam}->append_tell 
		: $aln;

	return $self->{_aln_idc};
}


=head2 remove_aln

Remove a Sam::Alignment from Sam::Seq (including position bins, if it exits).
 Returns the removed Sam::Alignment or undef. 

=cut

sub remove_aln{
	my ($self, $id);
	
	my $aln = delete $self->{_alns}{$id};
	defined $aln || return; 
	
	my $bin = $self->_bin($aln);
	my $idx = List::Util::first {$self->{_bin_alns}[$_] == $id} 1..@{$self->{_bin_alns}};
	defined $idx || return;
	
	splice(@{$self->{_bin_scores}[$bin]}, $idx, 1);
	splice(@{$self->{_bin_alns}[$bin]}, $idx, 1);
	return $aln;
}



=head2 consensus

Calculate and the consensus sequence from state matrix. Returns a Fastq::Seq
 object, with an additional entry C<< $seq->{cov} >>, which contains phred
 like ascii coded string, representing the per base coverages, offset 33.

=cut

sub consensus{
	my $self = shift;
	my @hcrs = @_;
	# init state matrix
	$self->{_state_matrix} = [map{[]}1..$self->len];
	$self->_state_matrix(@hcrs);# unless $self->_state_matrix();
	$self->_add_pre_calc_fq(@hcrs) if @hcrs;
	$self->_consensus;
	return $self->{con};
}


=head2 variants

=cut

sub variants{
	my ($self) = @_;
	$self->_state_matrix();# unless $self->_state_matrix();
	return $self->_variants;
}

=head2 coverage

Calculate and return LIST of accurate per base coverages.

=cut

sub coverage{
	my ($self) = @_;
	# calculate from _state_matrix, not the fastest way but accurate and
	#  already implemented :)
	$self->_state_matrix();# unless $self->_state_matrix();
	my @covs;
	foreach my $col(@{$self->{_state_matrix}}){
		if($col){
			my $s = 0;
			$s += $_ for (grep{$_}@$col);
			push @covs, $s;
		}else{
			push @covs, 0;
		}
	};
	return @covs;
}

=DEPRECATED freq2phred
#=head2 freq2phred
#
#Uses square root of the coverage probability of the most prominent
# state to compute phred like score,ranging from 0 to 40.
#
#  p(cov)  phred   Base call accuracy
#  0.0      0.00    0.00 %
#  0.1     12.65   94.57 %
#  0.2     17.88   98.37 %
#  0.4     25.29   99.70 %
#  0.6     30.98   99.92 %
#  0.8     35.77   99.97 %
#  1.0     40.00   99.99 %
#
#=cut
#
#sub freq2phred{
#	my ($self, $freq) = @_;
#	return 0  unless $freq;
#	# probability of state
#	my $p = $freq/$self->{max_coverage};
#	$p = 1 if $p > 1;
#	return int((sqrt($p) * 40) + .5); # phred
#}
#
#=head2 phred2freq
#
#Takes phred ranging from 0 to 40 and returns probable frequency of
# the state, used for phred calculation.
#
#=cut
#
#sub phred2freq{
#	my ($self, $phred) = @_;
#	return 0 unless $phred;
#	# phred = sqrt($p) * 40 => (phred /40)**29 = $p 
#	my $p = ($phred/40)**2;
#	return int(($p * $self->{max_coverage}) + 0.5);
#}
#
#
#=head2 char2phred
#
#
#=cut
#
#sub char2phred{
#	
#}
#
#=head2 phred2char
#
#
#
#=cut
#
#sub phred2char{
#	my ($self, $phred) = @_;
#	return chr($phred+$self->{phred_offset});
#}
=cut

=DEPRECATED approx_coverage
#=head2 approx_coverage
#
#Calculate and returns a LIST of coverage values, one value every $BinSize
# window. Values are estimated by the number of alns in each bin plus the 
# sum of alignments from previous bins, overlapping the current bin. Read
# length is assumed to be 100.
#
#This approximation is much faster than the exact per base coverage 
# calculation based on the state matrix. It is used to determine which 
# (high coverage) regions to mask before the second pass of shrimp.
#
#=cut
#
#sub approx_coverage{
#	my ($self) = @_;
#	
#	my @tmp;
#	my @covs;
#	my @cs = (0,0,0,0,0);
#	my $cum;
#	foreach my $bin(@{$self->{_bin_alns}}){
#		my $c = scalar @$bin;
#		push @tmp, $c;
#		push @cs, $c;
#		$cum += ($c - shift @cs);		
#		push @covs,$cum;
#	}
#	return  @covs;
#}


# DEPRECATED
#=head2 high_coverage_regions
#
#Determine high coverage regions (HCR) using C<approx_coverage>. HCRs have a
# minimum length of C<< int(READ_LENGTH / $BinSize) * 2 >> a minimum mean
# coverage of 70% of the coverage cutoff and no bin with a coverage below
# 20% of the coverage cutoff.
#
#=cut
#
##TODO: Values are educated guesses and only work for cov_co 50 and readlength 100
#
#sub high_coverage_regions{
#	my ($self) = @_;
#	my @covs = $self->approx_coverage();
#	
#	# debug -> take up the coverage
#	# @covs = map{$_*7}@covs;
#	
#	my @hcr;
#	my $hcrl = 0; # length
#	my $hcrs = 0; # sum
#	my $c;
#	
#	# run through coverage, search for HCRs
#	for(my $i=0; $i<@covs; $i++){
#		$c = $covs[$i];
#		print $c," ";
#		# hcr
#		if($c > 35 || ($c > 10 && (($hcrs+$c)/($hcrl+1)) > 0.7)){
#			$hcrs+=$c;
#			$hcrl++;
#			next;
#		}
#		# no_hcr
#		if($hcrl > 8){
#			# save the hcr, keep 4 bins at start and end unmasked for overlapping reads
#			push @hcr, [$i-$hcrl+4, $hcrl-4];
#		}
#		# reset hcr
#		$hcrl && ($hcrl = 0);
#		$hcrs && ($hcrs = 0);
#	}
#	
#	# end of seq, check once more for hcr
#	if($hcrl > 8){
#		# save the hcr, keep 4 bins at start unmasked, not at the end since
#		# obviously the last bins werent empty
#		push @hcr, [@covs-$hcrl+4, $hcrl];
#	}
#	print "\n";
#	return @hcr;
#}
=cut

=head2 chimera

Simple chimera detection. Searches for two hits of identical 
 reads in opposite orientation on the same reference sequences. Calculates
 a breakpoint for the chimera folding based on the distance between 
 multiple such twice mapping reads and returns the exstimated position
 of this breakpoint and the stddev from this point, in case enough evidence 
 has been found ($MaxCoverage/4 twice mapping reads) or undef otherwise.
 
           |
  0--->        10<---
     3<---   8--->


=cut

sub check_chimera{
	my $self = shift;
	my @alns = $self->alns();
	
	my %alns_by_id;
	my @breakpoints;
	my $bp_sum;
	
	# 
	foreach my $aln (@alns){
		push @{$alns_by_id{$aln->qname}}, $aln;
	}
	
	while(my ($k, $v) = each %alns_by_id){
		# two hits
		# three hits is strange...
		unless (@$v == 2){
			delete $alns_by_id{$k}; 
			next;		
		}; 

		# opposite direction
		my ($aln1, $aln2) = sort{$a->pos <=> $b->pos}@$v; # sort alns by pos
		unless ($aln1->is(REVERSE_COMPLEMENT) xor $aln2->is(REVERSE_COMPLEMENT)){
			delete $alns_by_id{$k};
			next;
		}
		my $bp = $aln1->pos + (($aln2->pos - $aln1->pos +length($aln1->seq))/2); 
		push @breakpoints, $bp;
		$bp_sum+=$bp;
	}
	
	return undef unless @breakpoints > $MaxCoverage/4;

	print Dumper(\@breakpoints);
	
	my $bp_mean = $bp_sum/@breakpoints;
	my $bp_dev = _stddev(\@breakpoints, $bp_mean);
	# round mean
	return (int($bp_mean + 0.5), int($bp_dev + 0.5)); 
}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 id

Get/Set the file handle.

=cut

sub id{
	my ($self, $id) = @_;
	$self->{id} = $id if $id;
	return $self->{id};
}

=head2 len

Get/Set the length.

=cut

sub len{
	my ($self, $length) = @_;
	$self->{len} = $length if $length;
	return $self->{len};
}

=head2 ref

Get/Set the reference sequence Fasta::Seq or Fastq::Seq object.

=cut

sub ref{
	my ($self, $ref) = @_;
	$self->{ref} = $ref if $ref;
	return $self->{ref};
}

=head2 con

Get the consensus sequence, a Fastq::Seq object with an additional entry 
 C<< $seq->{cov} >>, which contains an phred like ascii coded string, 
 representing the per base coverages, offset 33. Calls the consensus method
 in case no con object is present.

=cut

sub con{
	my ($self) = @_;
	$self->consensus unless $self->{con};
	return $self->{con};
}


=head2 next_aln

Returns the next Sam::Alignment of Sam::Seq. The alns are internally stored 
 in a hash and retrieved using each. Therefore the method behaves like each:
 It returns an empty list or undef, respectively, once after all alignments 
 have been returned. It returns alignments in apparently random order, 
 consistent as long as no alignments are added or removed. 

C< get_alns > returns in identical order as long as no alignments
 are added or removed.
 

=cut

sub next_aln{
	my ($self) = @_;
	my $aln;
	while(1){
		if($self->{sam}){
			my $pos = (each %{$self->{_alns}})[1];
			return undef unless defined $pos;
			$aln = $self->{sam}->aln_by_pos($pos);	#
			
		}else{
			$aln = (each %{$self->{_alns}})[1];
		}
		return undef unless defined($aln);
		last if &{$self->{_is}}($aln);
	}
	return $aln
}


=head2 alns(<SORTED_BY_POS>)

Returns number of Sam::Alignments of the Sam::Seq in scalar context, a list
 of all Sam::Alignments in list context. Order is identical to C<next_aln> as
 long as no alignments are added or removed or sorted by position, if first
 argument to call is TRUE.

  @alns = $pb->alns();     # LIST of alignments
  $num_alns = $pb->alns(); 
    # faster than getting list and using it in scalar context
  @alns_sorted_by_pos = $pb->alns(1);
    # LIST of alignments, sorted by pos, slower than without sorting


=cut

sub alns{
	my ($self, $sorted_by_pos) = @_;
	wantarray || return scalar keys %{$self->{_alns}};
	if($self->{sam}){
		# get indices from _aln and retrieve objects from parser
		if($sorted_by_pos){
			return sort{ $a->{'pos'} <=> $b->{'pos'} } map{ $self->{sam}->aln_by_pos($_) }values %{$self->{_alns}};
		}else{
			return map{ $self->{sam}->aln_by_pos($_ ) }values %{$self->{_alns}};
		}
	}else{
		# return objects from _aln
		if($sorted_by_pos){
			return sort{ $a->{'pos'} <=> $b->{'pos'} } values %{$self->{_alns}};
		}else{
			return values %{$self->{_alns}};
		}

	}
}


=head2 alns_by_bins

Returns list of bins, each containing list of Sam::Alignments, decendingly 
 ordered by their score.

=cut

sub alns_by_bins{
	my $self = shift;
	my @bins;
	foreach my $bin (@{$self->{_bin_alns}}){
		push @bins, [map{
			$self->{sam}
				? $self->{sam}->aln_by_pos( $self->{_alns}{$_} )
				: $self->{_alns}{$_};
		}@$bin];
	};
	return @bins;
}


=head2 bin

Compute the bin of given alignment. The bin is determined based on the 
 position of the center of the read (read_start + read_length/2) and the
 currently set bin size. Returns the bin (INT).

=cut

sub bin{
	my ($self,$aln) = @_;
	return int(( $aln->pos + ( length($aln->seq)/2 )) / $self->{bin_size})
}


=head2 is

Get/Set conditions that determine which alignments are returned by the 
 next methods. Takes either a reference to a list of property bitmasks 
 or a code reference to a customized test which returns 1 and 0 respectively.
 For details on bitmasks see L<Sam::Alignment>.
 
The test routine is executed with the parameters C<$parser_obj, $aln_obj> 
 and for C<next_pair()> additionally with C< $aln_obj2 >.

  # parser returning only BAD_QUALITY alns
  my $sp = Sam::Parser->new(
  	is => [Sam::Alignment->BAD_QUALITY]
  );
  
  # customized parser that only returns reads with a GC content > 70%.
  my $sp = Sam::Parser->new(
  	is => sub{
  	my ($self, $aln) = @_;
  	return ($aln->seq =~ tr/GC//) / length($aln->seq) > .7 ? 1 : 0;
  })
  
  
  
=cut

sub is{
	my ($self, $is) = @_;
	if($is){
		if(CORE::ref($is) eq 'ARRAY'){
			$self->{_is} = eval 'sub{$_[0]->is('.join(', ', @$is).')}';
		}elsif(CORE::ref($is) eq 'CODE'){
			$self->{_is} = $is;
		}else{
			die (((caller 0)[3])." neither ARRAY nor CODE reference given!\n");
		}	
	}
	return $self->{_is};
}

##------------------------------------------------------------------------##

=head1 Private Methods

=head2 _init_read_bins

Initialize bin structure required for storage of mapped read scores and ids.

=cut

sub _init_read_bins{
	my $self = shift;
	my $last_bin = int($self->len / $self->{bin_size});
	$self->{_bin_scores} = [map{[]}(0..$last_bin)];
	$self->{_bin_alns} = [map{[]}(0..$last_bin)];
	$self->{_bin_lengths} = [map{[]}(0..$last_bin)];
	$self->{_bin_bases} = [(0) x ($last_bin+1)];
}

=head2 _state_matrix

=cut

sub _state_matrix{
	my $self = shift;
	# state matrix
	my @S = @{$self->{_state_matrix}};
	# predefined states
	my %states = %{$self->{_states}};
	
	while(my $aln = $self->next_aln){
		
		###################
		### prepare aln ###
		# get read seq
		my $seq = $aln->seq;
		my $orig_seq_length = length($seq);
		
		next unless $orig_seq_length > 50;
		
		# get read cigar, eg 80M2D3M1IM4
		my @cigar = split(/(\d+)/,$aln->cigar);
		shift @cigar;
		
		$V->exit("Empty Cigar") unless @cigar;
		
		# reference position
		my $rpos = $aln->pos-1;

		if($Trim){
			##################
			### InDelTaboo ###
			# this also removes leading/trailing InDels regardless of InDelTaboo
			# trim head
			my $mc = 0;
			my $dc = 0;
			my $ic = 0;
			my $InDelTabooLength = int($orig_seq_length * $InDelTaboo + 0.5);
			for(my $i=0; $i<@cigar;$i+=2){
				if($cigar[$i+1] eq 'M'){
					if($mc + $ic + $cigar[$i] > $InDelTabooLength){
						if($i){# there was something before this match
							# only cut before this match
							# trim cigar
							splice(@cigar, 0, $i);
							# adjust rpos
							$rpos+= ($mc+$dc);
							# trim seq
							substr($seq, 0, $mc+$ic, '');
						}
						last;
					}
					$mc+=$cigar[$i];
				}elsif($cigar[$i+1] eq 'D'){
					$dc+= $cigar[$i];
				}elsif($cigar[$i+1] eq 'I'){
					$ic+= $cigar[$i];
				}else{
					$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
				}
			}
			
			# have to have kept at least 50 bps and 70% of original read length
			#  to consider read for state matrix
			next if length($seq) < 50  || (length($seq)/$orig_seq_length) < 0.7;
			
			# trim tail
			my $tail=0;	
			for(my $i=$#cigar-1; $i;$i-=2){
				if($cigar[$i+1] eq 'M'){
					$tail+=$cigar[$i];
					if($tail > $InDelTabooLength){
						if($i < $#cigar-1){# there is after this match
							# only cut before this match
							my $tail_cut = $tail-$cigar[$i]; 
							# trim cigar
							splice(@cigar, -($#cigar-($i+1)));
							# trim seq
							substr($seq, -$tail_cut, $tail_cut, '');
						}
						last;
					}
				}elsif($cigar[$i+1] eq 'D'){ # ignore leading deletions, but adjust rpos
	
				}elsif($cigar[$i+1] eq 'I'){
					$tail+=$cigar[$i];
				}else{
					$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
				}
			}
	
			# have to have kept at least 50 bps and 70% of original read length
			#  to consider read for state matrix
			next if length($seq) < 50  || (length($seq)/$orig_seq_length) < 0.7;
		}
		
=DEPRECATED leading/trailing long indel removal
		#		# remove leading/trailing I
		#		# there should be no l/t D since the reads are aligned semiglobal
		#		# leading
		#		if($cigar[1] eq 'I'){
		#			#use Data::Dumper; print Dumper("leading I", $aln);
		#			substr($seq,0,$cigar[0],''); # adjust read
		#			splice(@cigar,0,2); # adjust cigar
		#		}
		#		# trailing
		#		my $e = $#cigar;
		#		if($cigar[$e] eq 'I'){
		#			#use Data::Dumper; print Dumper("trailing I", $aln);
		#			substr($seq,-$cigar[$e-1],$cigar[$e-1],''); # adjust read
		#			splice(@cigar,$e-1,2); # adjust cigar
		#		}
		#		
		#		# detect long I/Ds within first 3 M of read
		#		if($cigar[1] eq 'M' && $cigar[0] < 4 && @cigar > 3){
		#			if($cigar[3] eq 'I' && $cigar[2] > 2){
		#				#use Data::Dumper; print Dumper("leading long I", $aln);
		#				
		#				# "long" (>= 3bp) I within first 3 matches -> discard read start
		#				$rpos+=$cigar[0]; # increase rpos of M
		#				substr($seq,0,$cigar[0]+$cigar[2],''); # remove I and M from read
		#				splice(@cigar,0,4); # adjust cigar
		#			}elsif($cigar[3] eq 'D' && $cigar[2] > 2){
		#				#use Data::Dumper; print Dumper("leading long D", $aln);
		#				
		#				# "long" (>= 3bp) D within first 3 M -> discard read start
		#				$rpos+=($cigar[0]+$cigar[2]); # increase rpos by number of M + D
		#				substr($seq,0,$cigar[0],''); # remove M from read
		#				splice(@cigar,0,4); # adjust cigar
		#			}
		#		}
=cut		
		
		
		
		#######################
		### cigar to states ###
		my @states;
		
		# cigar counter, increment by 2 to capture count and type of cigar (10,M) (3,I) (5,D) ...
		
		my $qpos = 0;
		for(my $i=0; $i<@cigar;$i+=2){
			if($cigar[$i+1] eq 'M'){
				push @states, split(//,substr($seq,$qpos,$cigar[$i]));
				$qpos += $cigar[$i];
			}elsif($cigar[$i+1] eq 'D'){
				push @states, ('-') x $cigar[$i];
			}elsif($cigar[$i+1] eq 'I'){
				if($i){
					# append to prev state
					$states[$#states] .= substr($seq,$qpos,$cigar[$i]);
				}else{
					$states[0] = substr($seq,$qpos,$cigar[$i]);		
				}
				$qpos += $cigar[$i];
			}else{
				$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
			}
		}
		
		
		########################
		### states to matrix ###
		
		foreach my $state (@states){
			if (length ($state) > 1 && ! exists $states{$state}){
				$states{$state} = scalar keys %states; 
			}		
			($S[$rpos][$states{$state}])++;  # match/gap states always exist
			$rpos++;
		}

	
=DEPRECATED simultaneous cigar parsing and state matrix counting
		#	my $state; # buffer last match, required if followed by insertion
		#		for(my $i=0; $i<@cigar;$i+=2){
		#			if($cigar[$i+1] eq 'M'){
		#				my @subseq = split(//,substr($seq,0,$cigar[$i],''));
		#				foreach $_ (@subseq){
		#					($S[$rpos][$states{$_}])++;  # match states always exist
		#					$rpos++;
		#				}
		#				$state = $subseq[$#subseq];
		#			}elsif($cigar[$i+1] eq 'D'){
		#				for(1..$cigar[$i]){
		#					($S[$rpos][4])++;  # $states{'-'} is always 4 
		#					$rpos++;
		#				}
		#				$state = '-';
		#			}elsif($cigar[$i+1] eq 'I'){
		#				#unless ($state){print STDERR $aln->pos," : ",$rpos,"\n"} 
		#				my $complex_state;
		#				if($state){
		#					$complex_state = $state.substr($seq,0,$cigar[$i],'');
		#					($S[$rpos-1][$states{$state}])--; #
		#				}else{
		#					$complex_state = substr($seq,0,$cigar[$i],'');
		#				}
		#				# replace by complex state, add state idx to %states if new
		#				if(exists ($states{$complex_state})){
		#					#TODO: insertion before first M
		#					next if ($rpos-1 < 0);
		#					$S[$rpos-1][$states{$complex_state}]++
		#				}else{
		#					next if ($rpos-1 < 0);
		#					$states{$complex_state} = scalar keys %states;
		#					$S[$rpos-1][$states{$complex_state}]++;
		#				}
		#				#($S[$rpos-1][exists ($states{$complex_state}) ? $states{$complex_state} : $states{$complex_state} = keys %states])++; 
		#				#$seq[$#seq].= 
		#			}else{
		#				$V->exit("Unknown Cigar '".$cigar[$i+1]."'");
		#			}
		#		}
=cut

	}
	
	
	# return state matrix
	$self->{_state_matrix} = \@S;
	$self->{_states} = \%states;
	
	return $self;
	
}

=head2 _add_pre_calc_fq

Add partial, already corrected sequences information (FASTQ) to the 
 state_matrix, to include them in consensus.

=cut

sub _add_pre_calc_fq{
	my $self = shift;
	my @S = @{$self->{_state_matrix}};
	foreach my $coords (@_){
		my ($seq) = $self->ref->slice_seq($coords);
		my @seq = split (//, $seq->seq);
		my @freqs = Sam::Seq->Phreds2Freqs($seq->phreds);
		
		use Data::Dumper;
		for(my $i=0; $i<length($seq->seq); $i++){
			# never add 0, if nothing more matches, a 0 count might be introduced
			next unless $freqs[$i];
			die Dumper($self->ref, $coords) unless defined $self->{_states}{$seq[$i]};
			($S[$i+$coords->[0]][$self->{_states}{$seq[$i]}])+=	$freqs[$i];
		}
	}
	$self->{_state_matrix} = \@S;
	return $self;
}


=head2 _consensus

=cut

sub _consensus{
	my $self = shift;
	my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
	my $seq = '';
#	my $qual = '';
#	my $covs = '';
	my @freqs;
	my $col_c = -1;
	foreach my $col (@{$self->{_state_matrix}}){
		$col_c++;
		#my $cov = @$col;
		# uncovered col
		unless (scalar @$col){
			$seq.= $self->{ref} ? substr($self->{ref}{seq}, $col_c, 1) : 'n';
#			$qual.=chr(0+$self->{phred_offset});
#			$covs.=chr(0+$self->{phred_offset});
			push @freqs, 0;
			next;
		}
		
		my $idx=undef;
		my $max_freq=0;
		my $i;
		my $cov;
		
		# majority vote
		for($i=0; $i<@$col; $i++){
			# get all defined states
			if (defined(my $freq = $col->[$i])){
				# add state frequency to coverage
				$cov+=$freq; 
				# check if more frequent than previous states
				if($freq > $max_freq){
					# exception: long prominent inserts (>3) are very ugly,
					# probable mapping artefacts caused by cheap gap costs 
					# compared to missmatch which might lead to long gaps
					# at read ends close to error rich regions. These long
					# gaps should not be considered.
					# Insert state has to have idx > 4 (not A,T,G,C or -)
					# for check to make sense
					next if($i > 4 && length $states_rev{$i} > 3);
										
					$max_freq = $freq;
					$idx = $i; 
				};
			} 
		};
		
		# check $max_freq, necessary due to long gap exception
		unless ($max_freq){
			$seq.= $self->{ref} ? substr($self->{ref}{seq}, $col_c, 1) : 'n';
#			$qual.=chr(0+$self->{phred_offset});
#			$covs.=chr(0+$self->{phred_offset});
			push @freqs, 0;
			next;
		}
		
		# insertion on reference
		next if $idx == 4; 
		
		# get most prominent state
		my $con = $states_rev{$idx};
		$seq.= $con;
		# entropy based quality
		#$qual.= $self->_phred_Hx($col) x length($con);
		# coverage based quality
#		$qual.= chr($self->freq2phred($max_freq) + $self->{phred_offset}) x length($con);
#		$covs.= chr($cov + $self->{phred_offset}) x length($con);
		push @freqs, ($max_freq) x length($con);
	}

	$self->{con} = Fastq::Seq->new(
		'@'.$self->{id},
		$seq,
		'+',
		Fastq::Seq->Phreds2Char( [Sam::Seq->Freqs2Phreds(@freqs)] , $self->{phred_offset} ),
		cov => Fastq::Seq->Phreds2Char([@freqs], $self->{phred_offset}),
		phred_offset => $self->{phred_offset}
	);
	
	return $self;
}



=head2 _variants

=cut

sub _variants{
	my $self = shift;
	my %p = (
		min_prob => 0.1,
		accuracy => 5,
		@_
	);
	
	
	#print Dumper($self->{_state_matrix});
	my @seq;
	my %states_rev = reverse %{$self->{_states}}; # works since values are also unique
	
	foreach my $col (@{$self->{_state_matrix}}){
		# cov
		unless($col){
			push @{$self->{covs}}, 0;
			push @{$self->{vars}}, ['?'];
			push @{$self->{freqs}},[''];
			push @{$self->{probs}},[''];
			next;
		}
		my $cov;
		my %vars;
		# variants
		for(my $i=0; $i<@$col; $i++){
			if (defined(my $v = $col->[$i])){
				$cov+= $v; 
				$vars{$states_rev{$i}} = $v;
			};
		};
		push @{$self->{covs}}, $cov;
		my @vars = sort{$vars{$b} <=> $vars{$a}}keys %vars;
		my @freqs = @vars{@vars};
		my @probs = map{sprintf("%0.".$p{accuracy}."f", $_/$cov)}@freqs;
		my $k = grep{$_>= $p{min_prob}}@probs;
		$k--;
		push @{$self->{vars}}, $k >= 0 ? [@vars[0..$k]] : ['?'];
		push @{$self->{freqs}}, $k >= 0 ? [@freqs[0..$k]] : [''];
		push @{$self->{probs}}, $k >= 0 ? [@probs[0..$k]] : [''];
	}
	# rel freq

	return $self;
}



sub _phred_Hx{
	my ($self, $col) = @_;
	my $total = 0;
	my @states = grep{$_}@$col;
	$total += $_ for @states;
	my @Px = map{$_/$total}@states;
	my $Hx;
	$Hx -= $_ for map{$_ * log $_}@Px;
	$Hx = 1 if $Hx > 1;
	$Hx = 0 if $Hx < 0;
	# prevent phreds >40 if actual coverage > $MaxCoverage
	# TODO: just a sloppy fix...
	# TODO !!!
	$total = $self->{max_coverage} if $total > $self->{max_coverage};
	return chr(
		int(
		((1-$Hx)/($self->{max_coverage}/$total))
		*30)							# phred range - phred minimum 
		+ 10							# phred minimum
		+ $self->{phred_offset});  		# phred offset
}



=head2 _stddev

Takes a reference to a list of values and returns mean and stddev of the
 list. Additionally takes the mean of the list as second parameter, in it
 has been calculated before to speed up computation.

=cut

sub _stddev{
	my($values, $mean1) = (@_);
	#Prevent division by 0 error in case you get junk data
	return undef unless scalar @$values;
	
	# calculate mean unless given
	unless(defined $mean1){
		# Step 1, find the mean of the numbers
		my $total1 = 0;
		$total1 += $_  for @$values;
		my $mean1 = $total1 / (scalar @$values);
	}
	
	
	# find the mean of the squares of the differences
	# between each number and the mean
	my $total2 = 0;
	$total2 += ($mean1-$_)**2 for @$values;
	my $mean2 = $total2 / (scalar @$values);
	
	# standard deviation is the square root of the
	# above mean
	return sqrt($mean2);
}





##------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



