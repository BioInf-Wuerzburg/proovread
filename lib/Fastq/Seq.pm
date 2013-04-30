package Fastq::Seq;

# $Id$

use warnings;
use strict;

use overload
	'""' => \&string,
	'.' => \&cat;
	
# preference libs in same folder over @INC
use lib '../';
use Verbose;

	
our $VERSION = '0.10';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;


##------------------------------------------------------------------------##

=head1 NAME 

Fastq::Seq.pm

=head1 DESCRIPTION

Class for handling FASTQ sequences.

=head1 SYNOPSIS

=cut

=head1 CHANGELOG

=head2 0.10

=over

=item [BugFix] C<< $fq->qual-window >> sets end of last window correctly if
 it coincides with last position of sequence.

=item [Feature] Added C<< $fq->compelment >> and C<< $fq->reverse_compelment >>

=item [BugFix] C<< $fq->substr_seq >> dies on no arguments

=item [Feature] C<< Fasta::Seq->Cat >> concatenates seq objects, "." and ".="
 are overloaded with this method

=item [Change] C<< $fq->qual_seq >> and C<< $fq->qual_lcs >> only return
 coordinates, not the qual sequence since substr_seq now takes the third
 argument as a replacement.

=item [Change] C<< $fq->substr_seq >> replaces C<< $fq->slice_seq >>. It 
 fully supports the syntax of perls C<substr> including replacements.
 It appends the id by  .COUNTER and adds the coordinates as attribute 
 SUBSTR:OFFSET,LENGTH to the description.

=back

=head2 0.09

=over

=item [Feature] C<< $fq->phred_transfrom >> to change the phred offset of
 the quality seq to the specified value.

=item [Change] Phred conversion done with C<unpack>, which is much faster.

=back

=head2 0.08

=over

=item [Change] Id and description parsing regexs carry C</o> flag.

=item [Change] << $fq->new() >> now provides cloning functionality, if called on
 object instead of Class. Using of << $fq->clone >> is DEPRECATED, but kept for 
 backward compability, simply calling << $fq->new >>.

=back

=over 12

=item 0.07 [Thomas Hackl 2012-11-6]

Bugfix. C<< $fq->substr_seq, $fq->substr_qual_lcs >> now add "_O:<OFFSET>_L:<LENGTH>"
 to the C<id>, not to the C<seq_head> of each newly created object.

=item 0.06 [Thomas Hackl 2012-10-31]

Bugfix. Faulty regex caused C<< $fq->desc >> to fail.

=item 0.05 [Thomas Hackl 2012-10-02]

C<< $fq->substr_seq, $fq->substr_qual_lcs >> now add "_O:<OFFSET>_L:<LENGTH>"
 to the C<seq_head> of each newly created object.

=item 0.04 [Thomas Hackl 2012-10-01]

Added C<< $fq->substr_seq, $fq->substr_qual_lcs >>. Allows to extract custom
 stretches and lcs stretches from sequence object and return them as new
 objects.

=item 0.03

Change new() behaviour. First 4 parameter are 4 lines of FASTQ, additional 
 parameters follow in key => value format. Makes it more straight forward 
 to parse lines from file to object.

=item 0.02

Added Picard naming scheme support.

=item 0.01

Initial Alignment module. Provides Constructor and generic accessor 
 methods.

=back

=cut

=head1 TODO

=over

=item Synopsis

=item Tests

=item Version

=item seq_decr qual_id qual_descr

=back

=cut


##------------------------------------------------------------------------##

=head1 Class Attributes

=cut

=head2 $V

Verbose messages are handled using the Verbose.pm module. To 
 customize verbose message behaviour, overwrite the attribute with
 another Verbose object created with the Verbose module.

=cut

our $V = Verbose->new();

our $Base_content_scans = {
	'N' => sub{	return $_[0] =~ tr/'N'// },
	'A' => sub{	return $_[0] =~ tr/'A'// },
	'T' => sub{	return $_[0] =~ tr/'T'// },
	'G' => sub{	return $_[0] =~ tr/'G'// },
	'C' => sub{	return $_[0] =~ tr/'C'// },
}; 

# escape bracketed character class special chars \,^,-,[ and ]
(our $Qual_lcs_range = 'YZ[\]^_`abcdefghij') =~ s/([\\\^\-\[\]])/\\$1/g;
our $Qual_lcs_min_length = 50;
our $Qual_lcs_regex = qr/([$Qual_lcs_range]{$Qual_lcs_min_length,})/o;

our $Qual_low_range = 'ABCDEFGHIJKLMNOPQRSTUVWX';
# doesn't make sense to have min length here
#our $Qual_low_min_length = 1;
our $Qual_low_regex = qr/([$Qual_low_range]+)/o;

our $Qual_window_size = 10;

our $Qual_window_min = 25;

# phred transform 
my $o33 = join("", map{chr($_+33)}(1..60));
my $o64 = join("", map{chr($_+64)}(1..60));
$o33 =~ s/([\\\^\-\[\]\/])/\\$1/g;
$o64 =~ s/([\\\^\-\[\]\/])/\\$1/g;
our %Tr;
$Tr{64} = eval 'sub{ $_[0] =~ tr/'.$o33.'/'.$o64.'/;}';
$Tr{33} = eval 'sub{ $_[0] =~ tr/'.$o64.'/'.$o33.'/;}';

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

=head2 Complement

=cut

sub Complement{	
	my ($class, $seq) = @_;
	$seq =~ tr/ATGCatgc/TACGtacg/;
	return $seq; 
} 

=head2 Reverse_complement

=cut

sub Reverse_complement{
	my ($class, $seq) = @_;
	$seq =~ tr/ATGCatgc/TACGtacg/;
	return reverse $seq;
}

=head2 Add_base_content_scan

=cut

sub Add_base_content_scan{
	my ($class, $patt) = @_;
	return unless defined $patt;
	$Base_content_scans->{$patt} = eval 'sub{ return $_[0] =~ tr/'.$patt.'//; }';
}

=head2 Qual_lcs_min_length [50]

Get/set $Fastq::Seq::Qual_lcs_min_length, which determines the minimum length
 for the qual_lcs() method

=cut

sub Qual_lcs_min_length{
	my ($class, $length) = @_;
	if(defined($length)){
		$Qual_lcs_min_length = $length;
		$Qual_lcs_regex = qr/([$Qual_lcs_range]{$Qual_lcs_min_length,})/o;	
	};
	return $Qual_lcs_min_length;
}

=head2 Qual_lcs_range ['YZ[\\]^_`abcdefghij'] / Qual_low_range [ABCDEFGHIJKLMNOPQRSTUVWX]

Get/set $Fastq::Seq::Qual_lcs_range and $Fastq::Seq::Qual_lcs_range, which 
 determine the character range for qual_lcs() and qual_low(), respectively. 
 Either provide a range of quality letters as one STRING or two phreds 
 values and the current phred offset.
 
  # high quality illumina lcs, phred 25-42 for offset 64
  Fastq::Seq::Qual_lcs_range('YZ[\]^_`abcdefghij') 
    # or
  Fastq::Seq::Qual_lcs_range(25, 42, 64)

  # low quality illumina, phred 25-42 for offset 64
  Fastq::Seq::Qual_low_range('ABCDEFGHIJKLMNOPQRSTUVWX')
    # or
  Fastq::Seq::Qual_low_range(1, 24, 64)   


=cut

sub Qual_lcs_range{
	my $class = shift;
	if(@_ == 1){
		$Qual_lcs_range = shift;
		# escape bracketed character class special chars \,^,-,[ and ]
		$Qual_lcs_range =~ s/([\\\^\-\[\]])/\\$1/g;
		$Qual_lcs_regex = qr/([$Qual_lcs_range]{$Qual_lcs_min_length,})/o;	
	}elsif(@_ == 3){
		my ($from, $to, $offset) = @_;
		$Qual_lcs_range = join('', map{chr($_ + $offset)}$from..$to);
		# escape bracketed character class special chars \,^,-,[ and ]
		$Qual_lcs_range =~ s/([\\\^\-\[\]])/\\$1/g;
		$Qual_lcs_regex = qr/([$Qual_lcs_range]{$Qual_lcs_min_length,})/o;	
	}elsif(@_){
		die sprintf("%s: %s",(caller 0)[3], "wrong number of parameter (1 or 3)");
	}
	return $Qual_lcs_range;
}

sub Qual_low_range{
	my $class = shift;
	if(@_ == 1){
		$Qual_low_range = shift;
		# escape bracketed character class special chars \,^,-,[ and ]
		$Qual_low_range =~ s/([\\\^\-\[\]])/\\$1/g;
		$Qual_low_regex = qr/([$Qual_low_range]+)/o;	
	}elsif(@_ == 3){
		my ($from, $to, $offset) = @_;
		$Qual_low_range = join('', map{chr($_ + $offset)}$from..$to);
		# escape bracketed character class special chars \,^,-,[ and ]
		$Qual_low_range =~ s/([\\\^\-\[\]])/\\$1/g;
		$Qual_low_regex = qr/([$Qual_low_range]+)/o;	
	}elsif(@_){
		die sprintf("%s: %s",(caller 0)[3], "wrong number of parameter (1 or 3)");
	}
	return $Qual_low_range;
}


=head2 Qual_lcs_regex [qr/([$Qual_lcs_range]{$Qual_lcs_min_length,})/o] /
 Qual_low_regex [qr/([$Qual_low_range]+)/o]

Get/set $Fastq::Seq::Qual_lcs_regex and $Fastq::Seq::Qual_low_regex, which 
 is used for lcs parsing and low quality parsing respectively. Set these 
 expressions with caution. Use Qual_lcs_range(), Qual_lcs_min_length() and 
 Qual_low_range instead.

=cut

sub Qual_lcs_regex{
	my ($class, $re) = @_;
	$Qual_lcs_regex = $re if $re;
	return $Qual_lcs_regex;
}

sub Qual_low_regex{
	my ($class, $re) = @_;
	$Qual_low_regex = $re if $re;
	return $Qual_low_regex;
}

=head2 Qual_window_min [25]

Get/Set min quality value for quality sliding window.

=cut

sub Qual_window_min{
	my ($class, $re) = @_;
	$Qual_window_min = $re if $re;
	return $Qual_window_min;
}

=head2 Qual_window_size [10]

Get/Set window size for quality sliding window.

=cut

sub Qual_window_size{
	my ($class, $re) = @_;
	$Qual_window_size = $re if $re;
	return $Qual_window_size;
}


=head2 Phreds2Char

Convert a LIST of phred values to a qual string.

=cut

sub Phreds2Char{
	my ($class, $p, $po) = @_;
	return pack("W*", map{$_+=$po}@$p);
}

=head2 cat

Concatenate Fastq::Seq object with each other. Can be used as Class method 
 as well as Object method.
Returns a new object. Keeps the id and other attributes from the first 
 provided object.
If the third parameter in Class syntax, the second one in object syntax is
 set to TRUE, the operands are swapped.

Fasta::Seq overloads "." with this method.

  # Class method
  $fab = Fastq::Seq->cat($fa, $fb);
  $fba = Fastq::Seq->cat($fa, $fb, 1); # swap order
  
  # Object method
  $fab = $fa->cat($fb);
  $fba = $fa->cat($fb, 1); # swap order
  
  # Overload
  $fa.= $fb; # append $fb and overwrite $fa

=cut

sub cat{
	my $class = shift unless ref $_[0]; # class usage
	my ($s1, $s2, $swap) = @_;
	
	my $re = $s1->new; # clone
	if($swap){
		$re->seq( $s2->seq.$re->seq );
		$re->qual( $s2->qual.$re->qual );
	}else{
		$re->seq( $re->seq.$s2->seq	);
		$re->qual( $re->qual.$s2->qual );
	}
	
	return $re;
}



##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

Create a new FASTQ seq object. Either provide "seq_head", "seq", 
 "qual_head", and "qual" as first four parameter or the four lines as one 
 STRING. Additional parameter can be specified in key => value format. 
 C<phred_offset> defaults to undef. 

When used on a Fastq::Seq object, it acts as cloning method, also 
 allowing for additional parameter in key => value format. Does not
 perform deep copy, not required for this object.

  # construct new object
  $fq = Fastq::Seq->new(<4_line_fastq_string>, phred_offset => 33);
  # or
  fq = Fastq::Seq->new(
  	<seq_head_string>,
  	<seq_string>,
  	<qual_head_string>,
    <qual_string>,
  	phred_offset => 33,
  );
  
  # clone object and set a new header at the same time
  $fq_clone = $fq->new(seq_head => 'NEW_HEADER');

=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	
	# object method -> clone + overwrite
	if($class = ref $proto){ 
		return bless ({%$proto, @_}, $class);
	}else{ # init empty obj
		$class = $proto;
		$self = {
				seq_head => '@',
				seq => '',
				qual_head => '+',
				qual => '',
				phred_offset => undef,
			};
	}
	
	# class method -> construct + overwrite
	if(@_){
		if(@_%2){ # create object from string
			my %self;
			@self{'seq_head','seq','qual_head','qual'} = split(/\n/, shift);
			$self = {
				%self,
				phred_offset => undef,
				@_	# overwrite defaults
			};
		}else{ # create object from array
			$self = {
				seq_head => $_[0],
				seq => $_[1],
				qual_head => $_[2],
				qual => $_[3],
				phred_offset => undef,
				@_[4..$#_]	# overwrite defaults
			};
			chomp(@$self{'seq_head','seq','qual_head', 'qual'});	# make sure all seqs loose there trailing "\n"
		}
	}
	
	return bless $self, $proto;
}




##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2 reverse_complement

Reverse complement the sequence. Returns modified seq object.

=cut

sub reverse_complement{
	my $self = shift;
	$self->{seq} = ref($self)->Reverse_complement($self->{seq});
	$self->{qual} = reverse $self->{qual};
	return $self;
}

=head2 reverse_complement

Complement the sequence. Returns modified seq object

=cut

sub complement{
	my $self = shift;
	$self->{seq} = ref($self)->Complement($self->{seq});
	return $self;
}

=head2 trim2qual_lcs

Trims nucleotide and quality sequence of object to qual_lcs, according to 
 Qual_lcs_range/Qual_lcs_min_length settings and returns the modified object
 or undef if no lcs was found.

=cut

# TODO: qual_lcs does not return qual_seq anymore!!!!!!

sub trim2qual_lcs{
	my ($self) = @_;
	my ($lcs) = $self->qual_lcs();
	if ($lcs){
		$self->{seq} = substr($self->{seq}, $lcs->[0], $lcs->[1]);
		$self->{qual} = $lcs->[2];
		return $self;
	}else{
		return undef;
	}
}



=head2 qual_lcs(<SORTED_BY_OCCURANCE>)

Calculate OFFSET (position) and LENGTH of contingious substrings of 
 predefined phred range (Qual_lcs_range()) and minimum length
 (Qual_lcs_min_length()).
Returns list of AARRAY tuples ([OFFSET,LENGTH], [OFFSET,LENGTH],...). 
Order is controlled by boolean parameter, if FALSE (default), LIST is sorted
 by LENGTH, longest first, else if it is TRUE, LIST is in order of occurance
 of the substrings.

=cut

sub qual_lcs{
	my ($self, $sorted_by_occurance) = (@_);
	my $qs = $self->{qual};
	my @re;
	while($qs =~ /$Qual_lcs_regex/g){
		push @re, [pos($qs) - length($1), length($1)];
	}	
	return $sorted_by_occurance ? sort{$b->[1] <=> $a->[1]}@re : @re;
}

=head2 qual_low

Calculate OFFSET (position) and LENGTH of low quality substrings of 
 predefined phred range (Qual_low_range()).
Returns list of ARRAY tuples ([OFFSET,LENGTH], [OFFSET,LENGTH],...), 
 in order of occurance.

=cut

sub qual_low{
	my $qs = $_[0]->{qual};
	my @re;
	while($qs =~ /$Qual_low_regex/g){
		push @re, [pos($qs) - length($1), length($1)];
	}	
	return @re;
}

=head2 mask_seq

Maskes sequence with "N"s. Takes list of ARRAY tuples ([OFFSET,LENGTH],
 [OFFSET,LENGTH],...) as provided by qual_lcs() or qual_low(). Returns 
 modified object.

=cut

sub mask_seq{
	my $self = shift;
	# mask seq substrings (offset, length) with "N"s.
	substr($self->{seq}, $_->[0], $_->[1], "N"x$_->[1]) for @_;
	return $self;
}

=head2 mask_qual_low

Maskes low quality regions of with "N"s in the nucleotide sequence. 
 Convenience function that is basically the same as 
 C<$fq->mask_seq($fq->qual_low())>. Low quality is therefore defined by 
 Qual_low_range()/Qual_low_regex.

=cut

sub mask_qual_low{
	my $self = shift;
	# mask seq substrings (offset, length) with "N"s.
	substr($self->{seq}, $_->[0], $_->[1], "N"x$_->[1]) for $self->qual_low();
	return $self;
}

=head2 substr_qual_lcs(<BOOL>)

Substr stretches of C< qual_lcs() > from the object and return resulting
 new partial objects. Set <BOOL> TRUE to get LIST sorted by 
 occurance instead of length of the lcs. Convenience function, basically the
 same as C<$fq->substr_seq($fq->qual_lcs())>

=cut

sub substr_qual_lcs{
	my ($self, $sorted_by_occ) = @_;
	my @lcs = $self->qual_lcs($sorted_by_occ);
	return @lcs ? $self->substr_seq(@lcs) : ();
}

sub slice_qual_lcs{
	warn 'use of "slice_qual_lcs" is deprecated, use "substr_qual_lcs" instead';
	return shift->substr_qual_lcs(@_);
}

=head2 substr_seq

Substr sequence. Takes the same parameter as perls C<substr>, either plain 
 or as a LIST of ARRAYREFS, each containing a set of parameter to allow for 
 multiple operations at once. 

In the first case the objects sequence is modified, the description appended
 by SUBSTR:<OFFSET>,<LENGTH> and the modified object is returned.

In the second case, the object itself is not modified, but a LIST of modified
 clones will be returned. Tie ids are appended by <.CLONECOUNTER>.
Methods like qual_low or qual_lcs provide these kind of LIST of ARRAYREFS.

  $fq->substr_seq(5);   	# nt 5 to end
  $fq->substr_seq(5,-5);    # nt 5 to fifth nt from the end
  $fq->substr_seq(-5)		# last 5 nts
  $fq->substr_seq(-100,50)  # 50 nts, starting 100 from end
  $fq->substr_seq(10,5,"AAAA", 'IIII') # replace 5 nts starting at pos 10 
    with 4 "A"s and replace corresponding qual values.
  ($fq1, $fq2) = $fq->substr_seq([5,100], [200,500]);
    # nts 5 to 100 and 200 to 500


=cut

sub substr_seq{
	my $self = shift;
	if(!ref $_[0] || @_ == 1){
		my $fq = $self->new; # clone
		
		my ($o, $l, $r, $q) = ref $_[0] ? @{$_[0]} : @_;
					
		die __PACKAGE__."::substr_seq: Not enougth arguments\n" 
			unless defined ($o);
			
		
		# replace
		if(defined $r){
			die __PACKAGE__."::substr_seq: seq and qual replacement need to be of same length \n$r\n$q\n" 
				unless defined ($q) && length ($r) eq length ($q);
			$fq->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
			substr($fq->{seq}, $o, $l, $r);
			substr($fq->{qual}, $o, $l, $q);
		}elsif(defined $l){
			$fq->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
			$fq->seq( substr($fq->{seq}, $o, $l) );
			$fq->qual( substr($fq->{qual}, $o, $l) );
		}else{
			$fq->desc_append(sprintf("SUBSTR:%d", $o));
			$fq->seq( substr($fq->{seq}, $o) );
			$fq->qual( substr($fq->{qual}, $o) );
		}
		return $fq;
	}else{
		my @new_fqs;
		my $clone_c = 0;
		
		foreach (@_){
			$clone_c++;
			my $fq = $self->new; # clone
			
			$fq->id($fq->id.".$clone_c");
	
			my ($o, $l, $r, $q) = @$_;
			
			die __PACKAGE__."::substr_seq: Not enougth arguments\n" 
				unless defined ($o);
			
			# replace
			if(defined $r){
				die __PACKAGE__."::substr_seq: seq and qual replacement need to be of same length \n$r\n$q\n" 
					unless defined ($q) && length ($r) eq length ($q);
				$fq->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
				substr($fq->{seq}, $o, $l, $r);
				substr($fq->{qual}, $o, $l, $q);
			}elsif(defined $l){
				$fq->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
				$fq->seq( substr($fq->{seq}, $o, $l) );
				$fq->qual( substr($fq->{qual}, $o, $l) );
			}else{
				$fq->desc_append(sprintf("SUBSTR:%d", $o));
				$fq->seq( substr($fq->{seq}, $o) );
				$fq->qual( substr($fq->{qual}, $o) );
			}
			push @new_fqs, $fq; 
		}
		return @new_fqs;
	}
}

sub slice_seq{
	warn 'use of "slice_seq" is deprecated, use "substr_seq" instead';
	return shift->substr_seq(@_);
}


=head2 phreds

Return the phred values of the quality string accorting to specified offset.

=cut

sub phreds{
	my $qs = $_[0]->{qual};
	my $po = $_[0]->{phred_offset};
	return map{$_-=$po}unpack("W*", $qs);
}



=head2 phred_transform

Transform phred 33 to 64 or 64 to 33.

=cut

sub phred_transform{
	my ($self) = @_;
	my $po = $self->{phred_offset};
	unless (defined $po && ($po == 33 || $po == 64)){
		die "Unknown or unsupported phred offset, cannot transform" 
	}
	$po = $po == 33 ? 64 : 33; # switch po 
	&{$Tr{$po}}($self->{qual});
	$self->{phred_offset} = $po;
	return $self;
}



=head2 base_content

=cut

sub base_content{
	my $self = shift;
	my $patt = shift;
	__PACKAGE__->Add_base_content_scan($patt) unless exists $Base_content_scans->{$patt};
	return &{$Base_content_scans->{$patt}}($self->seq)
}

sub qual_window{
	my ($self, $sorted_by_occurance) = (@_);
	
	use Data::Dumper;
	# sliding window
	# each window starts/ends above $min, values below $min are allowed as long
	# as the mean of the window values is higher than $min.
	# Cannot create overlapping windows!!
	# $min: 5
	# 663663636666363
	# 66366
	#   [6636] <- wont show
	#       63666636
	# seq starts with high window
	
	my @s = $self->phreds;
	my @re;
#	my @we;
#	my @ws;
	my $wsum = 0;
	my $wlen = undef;
	my $wsize = $Qual_window_size;
	my $wmin = $Qual_window_min * $wsize;


	# init first window
	$wsum += $_ for @s[0..$wsize-1];
	unless($wsum < $wmin){ # init window start
		#print "high init\n";
		$re[0][0] = 0;
		$wlen+=$wsize;
	}; 

	# init vars
	my $i=$wsize;
	my $j=0;
	
	#print " "x($wsize-1), $wsum,"\n";
	# loop
	for(;$i < @s; $i++){
		$wsum+=($s[$i]- $s[$j]); # new window sum
		#print " "x$i, $wsum,"\n";
		if($wsum < $wmin){# low
			if($wlen){# high ends
				my $x = $i-1;
				while($s[$x] < $Qual_window_min){ # go back to last pos over min
					$x--;
					$wlen--;
				}

				# make sure, window is still long enough to be accepted
				unless($wlen < $wsize){ 
					$re[$#re][1] = $wlen;
				}else{ # if not, remove start
					pop @re;
				}
				
				$wlen = 0;
				
				## no overlap ##
				# inti next window w/o overlap
				$i+=($wsize-1);
				$j+=($wsize-1);
				last if $j+$wsize > $#s; # make sure there is still enough sequence left
				$wsum = $s[$j+1];
				$wsum += $_ for @s[$j+2..($j+$wsize)];
			}
		}else{# high
			unless($wlen){
				if($s[$j+1] >= $Qual_window_min){ # high init only on over min pos
					# high init
					#print "high init\n";
					$re[@re][0] = $j+1;
					$wlen = $wsize;
				}
			}else{# high expand
				$wlen++
			}
		}
		$j++;
	}
	
	# set high end to last pos if seq ended high.
	# print $wlen , "\t", $wsum, "\t", $wmin,"\n";
	if($wlen){
		# print "yeah\n";
		$re[$#re][1] = scalar @s;
	}
	
	return $sorted_by_occurance ? sort{$b->[1] <=> $a->[1]}@re : @re;	

}

=head2

Append description by given attributes. Attributed are appended with
 whitespace delimiter.

=cut

sub desc_append{
	my ($self, @attr) = @_;
	my $cur_desc = $self->desc();
	$self->desc(
		$cur_desc
			?	join(" ", $cur_desc, @attr)
			:	join(" ", @attr)
	);
	return $self->desc;
}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 seq_head

Get/Set the seq_head.

=cut

sub seq_head{
	my ($self, $seq_head) = @_;
	if ($seq_head){
		# make sure seq_head starts with @
		$seq_head =~ /^@/ || substr($seq_head,0,0,'@');
		$self->{seq_head} = $seq_head;
	};
	return $self->{seq_head};
}

=head2 seq

Get/Set the seq.

=cut

sub seq{
	my ($self, $seq) = @_;
	if (defined $seq){
		$seq =~ tr/\n//d;
		$self->{seq} = $seq; 
	};
	return $self->{seq};
}

=head2 qual_head

Get/Set the qual_head line. 

=cut

sub qual_head{
	my ($self, $qual_head) = @_;
	if($qual_head){
		# make sure seq_head starts with @
		$qual_head =~ /^\+/ || substr($qual_head,0,0,'+');
		$self->{qual_head} = $qual_head 
	}
	return $self->{qual_head};
}

=head2 qual

Get/Set the seq.

=cut

sub qual{
	my ($self, $qual) = @_;
	if (defined $qual){
		$qual =~ tr/\n//d;
		$self->{qual} = $qual; 
	};
	return $self->{qual};
}

=head2 phred_offset

Get/Set the phred_offset.

=cut

sub phred_offset{
	my ($self, $offset) = @_;
	$self->{phred} = $offset if defined $offset;
	return $self->{phred_offset};
}


=head2 string

Get entire sequence as FASTQ string.

=cut

sub string{
	my ($self) = @_;
	return sprintf("%s\n%s\n%s\n%s\n", @$self{qw(seq_head seq qual_head qual)});
}

=head2 id

Get/set the reads id. Includes picard naming scheme like numbers at 
 the end (/1, /2). This is important especially when dealing with
 PacBio data because they use a trailing /<number> as part of the actual 
 read id.

=cut

sub id{
	my ($self, $id) = @_;
	if(defined $id){
		$id =~ s/^@//;
		$self->{seq_head} =~ s/[^@\s]+/$id/o;
	}else{
		($id) = $self->{seq_head} =~ /^@?(\S+)/o;
	}
	return $id;
}

=head2 id_no_picard

Get the reads id without the picard naming scheme like numbers at 
 the end (/1, /2).

=cut

sub id_no_picard{
	my ($self, $id) = @_;
	if(defined $id){
		$id =~ s/^@//;
		$self->{seq_head} =~ s/[^@\s]+?(\/\d+)?(\s|$)/$id$1$2/o;
	}else{
		($id) = $self->{seq_head} =~ /^@(\S+?)(?:\/\d+)?(?:\s|$)/o;  
	}
	return $id;
}

=head2 picard

Get the reads picard number (1,2..) from its id (/1, /2) or undef
 if no number.

=cut

sub picard{
	my ($self) = @_;
	# parse picard from head, cache picard
	my ($picard) = $self->{seq_head} =~ /^@(?:\S+?\/)(\d+)?(?:\s|$)/o;  
	return $picard;
}


=head2 desc

Get/set the reads description.

=cut

sub desc{
	my ($self, $desc) = @_;
	if(defined $desc){
		$self->{seq_head} = $desc 
			? sprintf('@%s %s', $self->id, $desc)
			: '@'.$self->id;
	}else{
		($desc) = $self->{seq_head} =~ /^\S+\s(.+)$/o;
	}
	return $desc || '';
}

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



