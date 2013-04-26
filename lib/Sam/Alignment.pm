package Sam::Alignment;

use warnings;
use strict;

# $Id$

use overload '""' => \&raw;

# preference libs in same folder over @INC
use lib '../';



our $VERSION = '0.07';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;


=head1 NAME 

Sam::Alignment.pm

=head1 DESCRIPTION

Class for handling sam alignments.

=cut

=head1 SYNOPSIS

  use Sam::Alignment;
  
  # directly build sam aln object
  $aln = Sam::Alignment->new();
  
  # get object from parser
  $sp = Sam::Parser->new(file => <SAMFILE>);
  $aln = $sp->next_aln;
  
  # get values
  $aln->qname     # query name
  $aln->raw       # entry raw string
  $aln->opt('XX') # value of optional field with tag 'XX'
  
  # test alns bit mask
  $aln->is_paired
  $aln->is(SAM::Alignment->PAIRED);
  $aln->is($aln->PAIRED);
  
  use Sam::Alingment ':flags'
  
  $sam_aln_obj->is(PAIRED);
  
  # true if read is paired and unmapped
  $aln->is(PAIRED, UNMAPPED);
  # true if reads is either duplicate or bad quality
  $aln->is(DUPLICATE & BAD_QUALITY);

=cut

=head1 CHANGELOG

=head 0.07

=over

=item [Feature] C<< $aln->raw >> overloads "";

=item [Change] C<split(/\s/, $sam, 12)> is replaced by C<split("\t", $sam, 12)>
 for better performance

=item [Feature] Field names are stored in @Sam::Alignment::_Fieldsnames

=back

=head2 0.06

=over

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=over

=item 0.05

C<< $aln->opt >> can now also be used to savely manipulate the optional
 fields.

=item 0.04

C<raw> line now always contains trailing newline.

=item 0.03

Added support for optional fields at the end of sam lines (TAG:TYPE:VALUE)

=item 0.02

Added 'raw' line to the alignment object attributes. 

=item 0.01

Initial Alignment module. Provides Constructor, generic accessor 
 and conditional C<is> filtering methods as well as exportable 
 "is_flags"

=back

=cut


=head1 TODO

=over

=item Tests

=back

=cut


=head1 Class ATTRIBUTES

=cut

our @_Fieldsnames = qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual opt);


=head1 Constructor METHOD

=head2 new

Create a sam alignment object. Takes either a sam entry as as string (one 
 line of a sam file) or a key => value representation of the sam fields
 C<qname flag rname pos mapq cigar rnext pnext tlen seq qual opt>.
 While the first eleven fields are regular, C<opt> contains a string of all 
 the optional fields added to the line.
 
Returns a sam alignment object. For more informations on the sam format see
 L<http://samtools.sourceforge.net/SAM1.pdf>.
 

=cut

sub new{
	my $class = shift;
	my $self;
	
	if(@_ == 1){ # input is string to split
		my $sam = $_[0];
		chomp($sam);
		my %sam;
		@sam{@Sam::Alignment::_Fieldsnames} = split("\t",$sam, 12); 
		$self = \%sam;
		$self->{raw} = $sam."\n";
	}else{ # input is key -> hash structure
		$self = {
			qname => undef,
			flag => undef,
			rname => undef,
			'pos' => undef,
			mapq => undef,
			cigar => undef,
			rnext => undef,
			pnext => undef,
			tlen => undef,
			seq => undef,
			qual => undef,
			opt => undef,
			@_,
			_opt => undef
		};
		if($self->{raw}){
			$self->{raw}.= "\n" unless $self->{raw} =~ /\n$/; # just to be safe, make it have one newline.
		}else{
			$self->{raw} = join("\t", @$self{qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual)});
			$self->{raw}.= "\t".$self->{opt} if $self->{opt};
			$self->{raw}.="\n";
		}
	}
	# overwrite defaults
	
	return bless $self, $class;
}

=head1 Object METHODS

=cut

=head2 is / is_<property>

Test generic alignment properties encoded in the sam bitmask flag. 
 Takes bitmasks, returns 1 or 0 accordingly.
 If you test multipe bitmasks, there are two different things you  might 
 want to see. a) all bitmasks have to match (AND linked) or b) at 
 least one bitmask has to match (OR linked). The achieve a) simply provide 
 the bitmasks as array to the method. For b) provide combined bitmasks 
 C<MASK_1or2 = (MASK_1 & MASK_2)>.
 
Named bitmask (property masks) for specific properties are provided for 
 better readability and consistency. The usage of these constants instead of
 actual bitmasks is recommended. The properties are read-only constants. 
 They can be exported individually by name or all at once with 
 C<use Sam::Alignment ':flags'>

  #property => bitmask value
  PAIRED => 0x1,
  MAPPED_BOTH => 0x2,
  UNMAPPED => 0x4,
  SECOND_READ_UNMAPPED => 0x8,
  REVERSE_COMPLEMENT => 0x10,
  SECOND_READ_REVERSE_COMPLEMENT => 0x20,
  FIRST => 0x40,
  SECOND => 0x80,
  SECONDARY_ALIGNMENT => 0x100,
  BAD_QUALITY => 0x200,
  DUPLICATE => 0x400
  
  # Property bitmasks   
  $sam_aln_obj -> is(SAM::Alignment->PAIRED);
  $sam_aln_obj -> is($sam_aln_obj->PAIRED);
  use Sam::Alingment ':flags'
  $sam_aln_obj -> is(PAIRED);

  # true if read is paired and unmapped
  is(PAIRED, UNMAPPED);
  # true if reads is either duplicate or bad quality
  is(DUPLICATE & BAD_QUALITY);
 
The C<is_<property>> methods are convenience functions, mainly equivalent to 
 C<is(<property>)>. The main advantage is that apart from 0 and 1 they 
 return undef in case the test itself does not make sence, that is if a 
 unpaired read is tested for paired read properties like 
 C<FIRST, SECOND, BOTH_MAPPED> or a non-first read is tested for 
 C<SECOND_READ_UNMAPPED, SECOND_READ_REVERSE_COMPLEMENT>.

  is_paired()
  is_mapped_both()
  is_first()
  is_second()
  is_unmapped()
  is_second_read_unmapped()
  is_reverse_complement()
  is_second_read_reverse_complement()
  is_secondary_alignment()
  is_bad_quality()
  is_duplicate()

=cut

our @flag_names = qw(
	PAIRED
	MAPPED_BOTH 
	UNMAPPED 
	SECOND_READ_UNMAPPED 
	REVERSE_COMPLEMENT 
	SECOND_READ_REVERSE_COMPLEMENT
	FIRST
	SECOND
	SECONDARY_ALIGNMENT
	BAD_QUALITY
	DUPLICATE
);

# declare property flags as constants
use constant {
	PAIRED => 0x1,
	MAPPED_BOTH => 0x2,
	UNMAPPED => 0x4,
	SECOND_READ_UNMAPPED => 0x8,
	REVERSE_COMPLEMENT => 0x10,
	SECOND_READ_REVERSE_COMPLEMENT => 0x20,
	FIRST => 0x40,
	SECOND => 0x80,
	SECONDARY_ALIGNMENT => 0x100,
	BAD_QUALITY => 0x200,
	DUPLICATE => 0x400
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

=head1 Accessor METHODS

=head2 qname

Get/Set the qname.

=cut

sub qname{
	my ($self, $qname) = @_;
	$self->{qname} = $qname if $qname;
	return $self->{qname};
}

=head2 flag

Get/Set the flag.

=cut

sub flag{
	my ($self, $flag) = @_;
	$self->{flag} = $flag if defined ($flag);
	return $self->{flag};
}

=head2 rname

Get/Set the rname.

=cut

sub rname{
	my ($self, $rname) = @_;
	$self->{rname} = $rname if $rname;
	return $self->{rname};
}

=head2 pos

Get/Set the pos.

=cut

sub pos{
	my ($self, $pos) = @_;
	$self->{pos} = $pos if $pos;
	return $self->{pos};
}

=head2 mapq

Get/Set the mapq.

=cut

sub mapq{
	my ($self, $mapq) = @_;
	$self->{mapq} = $mapq if $mapq;
	return $self->{mapq};
}

=head2 cigar

Get/Set the cigar.

=cut

sub cigar{
	my ($self, $cigar) = @_;
	$self->{cigar} = $cigar if $cigar;
	return $self->{cigar};
}

=head2 rnext

Get/Set the rnext.

=cut

sub rnext{
	my ($self, $rnext) = @_;
	$self->{rnext} = $rnext if $rnext;
	return $self->{rnext};
}

=head2 pnext

Get/Set the pnext.

=cut

sub pnext{
	my ($self, $pnext) = @_;
	$self->{pnext} = $pnext if $pnext;
	return $self->{pnext};
}

=head2 tlen

Get/Set the tlen.

=cut

sub tlen{
	my ($self, $tlen) = @_;
	$self->{tlen} = $tlen if $tlen;
	return $self->{tlen};
}

=head2 seq

Get/Set the seq.

=cut

sub seq{
	my ($self, $seq) = @_;
	$self->{seq} = $seq if $seq;
	return $self->{seq};
}

=head2 qual

Get/Set the seq.

=cut

sub qual{
	my ($self, $qual) = @_;
	$self->{qual} = $qual if $qual;
	return $self->{qual};
}

=head2 raw

Get/Set the raw line. 

NOTE: Following v0.04, inludes trailing newline.

=cut

sub raw{
	my ($self, $raw) = @_;
	if($raw){ # guarantee newline
		$self->{raw} = $raw =~ /\n$/ ? $raw : $raw."\n";
	}
	
	# {raw} is reset if any other value is changed and will be generated 
	# if required
	unless($self->{raw}){
		$self->{raw} = join("\t", @$self{qw(qname flag rname pos mapq cigar rnext pnext tlen seq qual)});
		$self->{raw}.= "\t".$self->{opt} if $self->{opt};
		$self->{raw}.="\n";
	}
	
	return $self->{raw};
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
		# reset raw, is regenerated with updated opt if required		
		$self->{raw} = undef;
		
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
			while($self->{opt} =~ /(\w\w):(\w):([^\t]+)/g){
				$self->{_opt}{$1} = [$2, $3];
			};
		}
	}
	return undef unless exists $self->{_opt}{$tag};
	return wantarray 
		? @{$self->{_opt}{$tag}} # type, value
		: $self->{_opt}{$tag}[1] # value only
}



=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



