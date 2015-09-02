package Sam::Alignment;

use warnings;
use strict;

use overload
    'bool' => sub{1},
    '""' => \&string;

use constant {
    GAP => '', # or '-'
};

our $VERSION = '1.1.5';

=head1 NAME

Sam::Alignment

=head1 DESCRIPTION

Class for handling sam alignments.

=cut

=head1 SYNOPSIS

  use Sam::Alignment;

  # usually get object from parser
  $sp = Sam::Parser->new(file => <SAM>);
  $aln = $sp->next_aln;

  # get values
  $aln->qname     # query name
  $aln->raw       # entry raw string
  $aln->opt('XX') # value of optional field with tag 'XX'

  use Sam::Alingment ':flags'

  # true if read is paired and unmapped
  $aln->is(PAIRED, UNMAPPED);
  # true if reads is either duplicate or bad quality
  $aln->is(PCR_DUPLICATE & VENDOR_FAIL);

=head1 Class ATTRIBUTES

=head2 $InvertScores;

=cut

our $InvertScores;

=head1 Class METHODS

=cut

=head2 InvertScores

Get/Set $Sam::Seq::InvertScore. Default OFF.

=cut

sub InvertScores{
	my ($class, $flag) = @_;
        return $InvertScores unless defined $flag;
	$InvertScores = $flag ? 1 : 0;
	return $InvertScores;
}


=head1 Constructor METHOD

=head2 new

Create a sam alignment object. Takes either a sam entry as as string (one
 line of a sam file) or a key => value representation of the sam fields
 C<qname flag rname pos mapq cigar rnext pnext tlen seq qual opt>.

Returns a sam alignment object. For more informations on the sam format see
 L<http://samtools.sourceforge.net/SAM1.pdf>.

=cut

my @ATTR_SCALAR = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual opt);
my @ATTR_SCALAR_def = (qw(* 0 * 0 255 * * 0 0 * *), undef);

my %SELF;
@SELF{@ATTR_SCALAR} = @ATTR_SCALAR_def;

sub new{
	my $class = shift;
	my $self;

	if(@_ == 1){ # input is string to split
		my $sam = $_[0];
		chomp($sam);
		my %sam;
		@sam{@ATTR_SCALAR} = split("\t",$sam, 12);
		$self = \%sam;
                $self->{qual} //= "*";

	}else{ # input is key -> hash structure
            $self = {
                %SELF,
                @_,
                _opt => undef
            };
	}
	# overwrite defaults


	return bless $self, $class;
}


=head1 Object METHODS

=head2 qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual

Get/Set ...

  my $id = $aln->qname();
  $aln->seq("AATTATA");

=cut

# called at eof
sub _init_accessors{
    no strict 'refs';

    # generate accessors for cache affecting attributes
    foreach my $attr ( qw(pos cigar tlen seq) ) {
        next if $_[0]->can($attr); # don't overwrite explicitly created subs
        *{__PACKAGE__ . "::$attr"} = sub {
            if (@_ == 2){
                $_[0]->_reset_cached_values();
                $_[0]->{$attr} = $_[1];
            }
            return $_[0]->{$attr};
        }
    }

    # generate simple accessors closure style
    foreach my $attr ( @ATTR_SCALAR ) {
        next if $_[0]->can($attr); # don't overwrite explicitly created subs
        *{__PACKAGE__ . "::$attr"} = sub {
            $_[0]->{$attr} = $_[1] if @_ == 2;
            return $_[0]->{$attr};
        }
    }
}


=head2 string

Get stringified alignment.

=cut

sub string{
    my ($self) = @_;
    my $s = join("\t", @$self{@ATTR_SCALAR[0..$#ATTR_SCALAR-1]});
    $s.= "\t".$self->{opt} if $self->{opt};
    return $s."\n";
}

=head2 is

Test flag encoded alignment properties. Takes bitmasks, returns 1 or 0
 accordingly. For OR testing of multiple bitmasks, provide maskes as array. For
 AND testing of bitmasks combine bitmasks C<(MASK_1 & MASK_2)>. Export named
 bitmasks with C<use Sam::Alignment ':flags'>

  # Property bitmasks
  use Sam::Alingment ':flags'

  # true if read is paired and unmapped
  is(PAIRED, UNMAPPED);
  # true if reads is either duplicate or bad quality
  is(DUPLICATE & BAD_QUALITY);

  #property => bitmask value
  PAIRED =>          0x1
  PROPER_PAIR =>     0x2
  UNMAP =>           0x4
  MUNMAP =>          0x8
  REVERSE =>        0x10
  MREVERSE =>       0x20
  READ1 =>          0x40
  READ2 =>          0x80
  SECONDARY =>     0x100
  QCFAIL =>        0x200
  DUP =>           0x400
  SUPPLEMENTARY => 0x800

=cut

# deprecated
# The C<is_<property>> methods are convenience functions, mainly equivalent to
#  C<is(<property>)>. The main advantage is that apart from 0 and 1 they
#  return undef in case the test itself does not make sence, that is if a
#  unpaired read is tested for paired read properties like
#  C<FIRST, SECOND, BOTH_MAPPED> or a non-first read is tested for
#  C<SECOND_READ_UNMAPPED, SECOND_READ_REVERSE_COMPLEMENT>.
#
#   is_paired()
#   is_mapped_both()
#   is_first()
#   is_second()
#   is_unmapped()
#   is_second_read_unmapped()
#   is_reverse_complement()
#   is_second_read_reverse_complement()
#   is_secondary_alignment()
#   is_bad_quality()
#   is_duplicate()


our @flag_names = qw(
	PAIRED
	PROPER_PAIR           MAPPED_BOTH
        UNMAP                 UNMAPPED
	MUNMAP                SECOND_READ_UNMAPPED
	REVERSE               REVERSE_COMPLEMENT
        MREVERSE              SECOND_READ_REVERSE_COMPLEMENT
	READ1                 FIRST
	READ2                 SECOND
	SECONDARY             SECONDARY_ALIGNMENT
	QCFAIL                BAD_QUALITY
	DUP                   DUPLICATE
        SUPPLEMENTARY
);

# declare property flags as constants
use constant {
    PAIRED => 0x1,
    PROPER_PAIR => 0x2,     MAPPED_BOTH => 0x2,
    UNMAP => 0x4,           UNMAPPED => 0x4,
    MUNMAP => 0x8,          SECOND_READ_UNMAPPED => 0x8,
    REVERSE => 0x10,        REVERSE_COMPLEMENT => 0x10,
    MREVERSE => 0x20,       SECOND_READ_REVERSE_COMPLEMENT => 0x20,
    READ1 => 0x40,          FIRST => 0x40,
    READ2 => 0x80,          SECOND => 0x80,
    SECONDARY => 0x100,     SECONDARY_ALIGNMENT => 0x100,
    QCFAIL => 0x200,        BAD_QUALITY => 0x200,
    DUP => 0x400,           DUPLICATE => 0x400,
    SUPPLEMENTARY => 0x800,
    # drives ncscoring penality for short alns, see _ncscore for details
    NCSCORE_CONSTANT => 40,
};


# export options
use Exporter qw(import);
our @EXPORT_OK = (@flag_names, '@flag_names');
our %EXPORT_TAGS = ('flags' => [@flag_names]);

# generic is() function
sub is{
	my $self = shift;
	foreach (@_){
		return 0 unless ($self->{flag} & $_);
	};
	return 1;
}

sub is_paired{
	return $_[0]->is(PAIRED);
}

sub is_mapped_both{
	my $self = shift;
	return undef unless $self->is(PAIRED);
	return $self->is(MAPPED_BOTH);
}

sub is_first{
	my $self = shift;
	return undef unless $self->is(PAIRED);
	return $self->is(FIRST);
}

sub is_second{
	my $self = shift;
	return undef unless $self->is(PAIRED);
	return $self->is(SECOND);
}

sub is_unmapped{
	my $self = shift;
	$self->is(UNMAPPED);
}

sub is_second_read_unmapped{
	my $self = shift;
	return undef unless $self->is(PAIRED, FIRST);
	return $self->is(SECOND_READ_UNMAPPED);
}

sub is_reverse_complement{
	my $self = shift;
	return $self->is(REVERSE_COMPLEMENT);
}

sub is_second_read_reverse_complement{
	my $self = shift;
	return undef unless $self->is(PAIRED, FIRST);
	return $self->is(SECOND_READ_REVERSE_COMPLEMENT);
}

sub is_secondary_alignment{
	my $self = shift;
	return $self->is(SECONDARY_ALIGNMENT);
}

sub is_bad_quality{
	my $self = shift;
	return $self->is(BAD_QUALITY);
}

sub is_duplicate{
	my $self = shift;
	return $self->is(DUPLICATE);
}


=head2 opt

Get/Set optional fields of TAG:TYPE:VALUE format in the sam file.
 Takes a either a single TAG id of an optional field as input and returns
 in SCALAR context the respective field VALUE, in ARRAY context TYPE and
 VALUE of the field.

 Without parameter, the string containing all optional fields is returned.

 Provide TAG, TYPE, VALUE to set an optional field. For possible TYPES see
 SAM specification.

 Provide TAG, TYPE, undef to remove an optional field.


=cut

sub opt{
	my ($self, $tag, $type, $value) = @_;
	# set
	if($type){

		# make sure, opt has been parsed, so it can be overwritten
		unless ($self->{_opt}){
			while($self->{opt} =~ /(\w\w):(\w):([^\t]+)/g){
				$self->{_opt}{$1} = [$2, $3];
			};
		}

		# remove
		unless(defined $value){
			delete $self->{_opt}{$tag};
		# add
		}else{
			$self->{_opt}{$tag} = [$type, $value];
		};

		my @opt;

		foreach my $tag (sort keys %{$self->{_opt}}	){
			push @opt , $tag.':'.$self->{_opt}{$tag}[0].':'.$self->{_opt}{$tag}[1];
		}

		$self->{opt} = join("\t", @opt);
	# get
	}else{
		return $self->{opt} unless $tag;
		return undef unless $self->{opt};
		unless ($self->{_opt}){
                    foreach(split("\t", $self->{opt})){
                        $self->{_opt}{substr($_, 0, 2)} = [substr($_, 3, 1), substr($_, 5)];
                    };
		}
	}
	return undef unless exists $self->{_opt}{$tag};
	return wantarray
		? @{$self->{_opt}{$tag}} # type, value
		: $self->{_opt}{$tag}[1] # value only
}


=head2 full_length

Get the length of the full query sequence including clipped bases. This will
crash if mapping software fails to annotate clipped alignments properly with
H/S.

=cut

sub full_length{
    my ($self) = @_;
    return $self->{full_length} if defined($self->{full_length});

    my $cigar = $self->cigar;
    my $l;
    while ($cigar =~ /(\d+)[MDHS]/g) { $l += $1; };

    $self->{full_length} = $l;
    return $l;
}


=head2 length

Get the length of the actually aligned sequence, meaning the number of query
bases in the alignment. For global alignments or hard clipping, this is equal to
length($aln->seq).

This can be a little bit slow as cigar string parsing might be
required. However, the result is cached, making successive calls less expensive.

=cut

sub length{
    my ($self) = @_;
    return $self->{length} if defined($self->{length});

    my $cigar = $self->cigar;
    my $l;
    if ($self->seq eq "*" || $cigar =~ /^\d+S|S$/) { # seq NA or soft clipped
        while ($cigar =~ /(\d+)[MD]/g) { $l += $1; };
    }else {
        $l = length($self->seq);
    }

    $self->{length} = $l;
    return $l;
}


=head2 seq_aligned

Get the aligned sequence, with D as gaps and I/S removed from seq. This allows
to correctly extract bases using reference positions.

=cut

sub seq_aligned{
    my ($self) = @_;

    return undef if $self->seq eq "*";

    my $cigar = $self->cigar;
    my $seq = $self->seq;
    my $aseq = '';
    my $pos = 0;
    if ($cigar =~ /(\d+)S/) {
        $pos+=$1; # account for softclip
    }

    while ($cigar =~ /(\d+)([MDIX=])/g) {
        if ($2 eq "I") {
            $pos+= $1
        } elsif ($2 eq "D") {
            $aseq.= "-" x $1;
        }  else {
            $aseq.= substr($seq, $pos, $1); $pos+=$1
        }
    };
    return $aseq;
}

=head2 seq_states

Get the states of a sequence based on cigar. By default gaps are represented by
char defined with GAP constant (''), can be overwritten with (gap => '-');

=cut

sub seq_states{
    my $self = shift;

    return undef if $self->seq eq "*";

    my %p = (gap => GAP, @_);

    my $cigar = $self->cigar;
    my $seq = $self->seq;
    my @s = ();
    my $pos = 0;
    if ($cigar =~ /^(\d+)S/) {
        $pos+=$1; # account for softclip
    }

    while ($cigar =~ /(\d+)([MDIX=])/g) {
        if ($2 eq "I") {
            $s[-1].= substr($seq, $pos, $1) if @s;
            $pos+=$1;
        } elsif ($2 eq "D") {
            push @s, ($p{gap}) x $1;
        } else {
            push @s, split(//, substr($seq, $pos, $1));
            $pos+=$1;
        }
    };
    return @s;
}

=head2 score, nscore, ncscore

  score:    AS:i
  nscore:   score/length
  ncscore:  score/length * CF

If $InvertScore is true, scores are multipled by -1. That way one can work
with/filter scores where originally smaller values are better, e.g. blasr
scores.

CF is a correction factor accounting for increased uncertainty in short
alignments. ncscore asymptotically approaches nscore for long alignments, but
penalizes shorter ones, based on NSCORE_CONSTANT [40].

  cf = length / ( NSCORE_CONTANT + length )

  |    bp | cf[10] | cf[20] | cf[40] |
  |-------+--------+--------+--------|
  |    10 |   0.50 |   0.33 |   0.20 |
  |    50 |   0.83 |   0.71 |   0.56 |
  |   100 |   0.91 |   0.83 |   0.71 |
  |   200 |   0.95 |   0.91 |   0.83 |
  |   500 |   0.98 |   0.96 |   0.93 |
  |  1000 |   0.99 |   0.98 |   0.96 |
  |  5000 |   1.00 |   1.00 |   0.99 |
  | 10000 |   1.00 |   1.00 |   1.00 |
  | 50000 |   1.00 |   1.00 |   1.00 |

=cut

sub score{
    my ($self) = @_;
    my $score = $self->opt("AS");
    return unless defined $score;
    return $InvertScores ? $score * -1 : $score;
}

sub nscore{
    my ($self) = @_;
    my $score = $self->score;
    return unless defined $score;

    return $score / $self->length;
}

sub ncscore{
    my ($self) = @_;
    my $nscore = $self->nscore;
    return unless defined $nscore;

    return $nscore * ($self->length/(NCSCORE_CONSTANT + $self->length));
}

##----------------------------------------------------------------------------##

=head1 Privates METHODS

=head2 _reset_cached_values

=cut

sub _reset_cached_values{
    $_[0]->{length} = undef;
    $_[0]->{full_length} = undef;
}

=head1 Aliases

=head2 raw = string

For backward compatibility or lazyness.

=cut

{ # alias for backward comp.
    no warnings 'once';
    *raw = \&string;
}

# init auto-accessors at eof to prevent any overwrites
__PACKAGE__->_init_accessors();


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut

1;
