package Fasta::Seq;

use warnings;
use strict;

# $Id$

# preference libs in same folder over @INC
use lib '../';

use Verbose;

use overload
	'.' => \&cat,
	'""' => \&string;


our $VERSION = '0.06';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;

##------------------------------------------------------------------------##

=head1 NAME 

Fasta::Seq.pm

=head1 DESCRIPTION

Class for handling FASTA sequences.

=head1 SYNOPSIS

=cut

=head1 CHANGELOG

=head2 0.06

=over

=item [BugFix] Regex to parse FASTA record doesn't fail anymore on 
 entries without description.

=item [BugFix] C<< $fq->substr_seq >> dies on no arguments

=item [Change] Refactored C<new()>. Handles creation of "empty" objects and 
 has some adjustments to head/id/desc parsing/generation.

=item [Change] Removed C<< $fa->clone >>. Functionality is provided by
 C<< $fa->new >>.

=item [Feature] C<< Fasta::Seq->Cat >> concatenates seq objects, "." and ".="
 are overloaded with this method

=item [Change] C<< $fq->substr_seq >> replaces C<< $fq->slice_seq >>. It 
 fully supports the syntax of perls C<substr> including replacements.
 It appends the id by  .COUNTER and adds the coordinates as attribute 
 SUBSTR:OFFSET,LENGTH to the description and allows multiple operations at
 once provided a LIST of ARRAYREFs containing the coordinates.

=back

=item [Feature] Added "cloning" behaviour to C<< $fa->new >>.

=item [Feature] Added class methods C<< Fastq::Seq->Reverse_complement >> 
 and C<< Fastq::Seq->Complement >>, as well as object methods 
 C<< $fa->reverse_complement >> and C<< $fa->complement >>.

=back

=head2 0.05

=over

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=over 12

=item 0.04 [Thomas Hackl 2012-10-31]

Change C<< Fastq::Seq->new >> parameter. Takes still a single string
 or additionally now a complete key => value construct.

=item 0.03 [Thomas Hackl 2012-10-28]

POD correction, C<< Fasta::Seq->new >> Synopsis.

=item 0.02 [Thomas Hackl 2012-10-24]

Added C<< byte_offset >> attribute to Fasta::Seq.

=item 0.01

Initial module. Provides Constructor and generic accessor 
 methods.

=back

=cut


=head1 TODO

=over

=item Synopsis

=item Tests

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

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut

=head2 Add_base_content_scan

=cut

sub Add_base_content_scan{
	my ($class, $patt) = @_;
	return unless defined $patt;
	$Base_content_scans->{$patt} = eval 'sub{ return $_[0] =~ tr/'.$patt.'//; }';
}

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

##------------------------------------------------------------------------##


=head1 Constructor METHOD

=head2 new

Create a new FASTA seq object. Either provide a single STRING 
 containing one FASTA record or a key => value construct, for which
 either C<seq_head> or C<id> is required.

  my $seq = ">seq1 blub\nATGC\n";
  Fasta::Seq->new($seq);
  # or
  Fasta::Seq->new(
  	id => 'seq1',
  	desc => 'blub',
  	seq => "ATGC",
  );
  # or
  Fasta::Seq->new(
  	seq_head => '>seq1 blub'
  	seq => "ATGC"
  );


=cut

sub new{
	my $proto = shift;
	my $self;
	my $class;
	# object method -> clone + overwrite
	if($class = ref $proto){
		$self = bless ({%$proto, @_}, $class)
	}else{ # init emtpy obj
		$class = $proto;
		$self = {
				seq_head => '',
				id => '',
				desc => '',
				seq => '',
				byte_offset => undef,
			};
	}
	
	if(@_){
		if(@_%2){ # input is string to split
			my %self;
			@self{'id','desc', 'seq'} = shift =~ m/
				(?:>?(\S*))			# id, >? for records
				(?:[^\S\n]([^\n]+))?\n	# desc, optional
				(.+)				# seq
			/xs;					# . includes \n
			$self = {
				%$self,
				%self,
				@_	# overwrite defaults
			};
			$self->{seq_head} = '>'.$self->{id};
			if($self->{desc}){
				$self->{seq_head}.= ' '.$self->{desc}
			}else{
				$self->{desc} = '';
			}
	
		}else{
			$self = {
				%$self,
				@_	# overwrite defaults
			};
			
			# make sure, id has not leading '>'
			$self->{id} =~ s/^>// if $self->{id};
			# make sure there is a seq_head entry
			if(!$self->{seq_head} && $self->{id}){
				$self->{seq_head} = '>'.$self->{id};
				$self->{seq_head}.= ' '.$self->{desc} if $self->{desc};
			}elsif($self->{seq_head}){
				my($id,$desc) = $self->{seq_head} =~ m/
					(?:>?(\S*))			# id, >? for records
					(?:[^\S\n]([^\n]+))?		# desc, optional
				/xs;	
				$self->{id} = $id; 
				$self->{desc} = $desc || '';
				
			}
		}
		# make sure, head has leading '>'
		substr($self->{seq_head},0,0,'>') unless $self->{seq_head} =~ /^>/;
		chomp($self->{seq});
		chomp($self->{seq_head});
		$self->{seq} =~ tr/\n//d; 	# remove all newlines from seq
	}
	
#	unless($self->{id} || $self->{seq}){
#		warn "Creation of incomplete FASTA entry: ".
#			($self->{id} ? "sequence missing" : "id missing");
#	}

	return bless $self, $class;
}







##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2 reverse_complement

Reverse complement the sequence.

=cut

sub reverse_complement{
	my $self = shift;
	$self->{seq} = ref($self)->Reverse_complement($self->{seq});
}

=head2 reverse_complement

Complement the sequence.

=cut

sub complement{
	my $self = shift;
	$self->{seq} = ref($self)->Complement($self->{seq});
}


=head2 cat

Concatenate Fasta::Seq object with either another Fasta::Seq object or a 
 plain STRING (sequence only). Can be used as Class method as well as Object 
methods.
Returns a new object. Keeps the id and other attributes from the first 
 provided object.
If the third parameter in Class syntax, the second one in object syntax is
 set to TRUE, the operands are swapped.

Fasta::Seq overloads "." with this method.

  # Class method
  $fab = Fasta::Seq->cat($fa, $fb);
  $fba = Fasta::Seq->cat($fa, $fb, 1); # swap order
  $faATGC = Fasta::Seq->cat($fa, 'ATGC'); # append plain sequence
  $fATGCa = Fasta::Seq->cat('ATGC', $fa); # prepend plain sequence
  
  # Object method
  $fab = $fa->cat($fb);
  $fba = $fa->cat($fb, 1); # swap order
  $fba = $fa->cat('ATGC', 1); # append + swap -> prepend plain sequence
  
  # Overload
  $fATGCa = 'ATGC'.$fb; # prepend
  $fa.= $fb; # append by $fb and overwrite $fa

=cut

sub cat{
	my $class = shift unless ref $_[0]; # class usage
	my ($s1, $s2, $swap) = @_;
	unless(ref $s1){
		die 'At least one operand has to be as Fasta::Seq object' unless ref $s2;
		($s1,$s2) = ($s2,$s1);
		$swap = $swap ? 0 : 1; # toggle swap
	}
	
	my $re = $s1->new; # clone
	if($swap){
		$re->seq(ref $s2
			? $s2->seq.$re->seq
			: $s2.$re->seq
		);
	}else{
		$re->seq(ref $s2
			? $re->seq.$s2->seq
			: $re->seq.$s2
		);
	}
	
	return $re;
}


=head2 base_content

=cut

sub base_content{
	my $self = shift;
	my $patt = shift;
	__PACKAGE__->Add_base_content_scan($patt) unless exists $Base_content_scans->{$patt};
	return &{$Base_content_scans->{$patt}}($self->seq)
}


=head2 substr_seq

Substr sequence. Takes the same parameter as perls C<substr>, either plain 
 or as a LIST of ARRAYREFS, each containing a set of parameter to allow for 
 multiple operations at once. 

In the first case the objects sequence is modified, the description appended
 by SUBSTR:<OFFSET>,<LENGTH> and the modified object is returned.

In the second case, the object itself is not modified, but a LIST of modified
 clones will be returned. The ids are appended by <.CLONECOUNTER>.
Methods like qual_low or qual_lcs provide these kind of LIST of ARRAYREFS.

  $fq->substr_seq(5);   	# nt 5 to end
  $fq->substr_seq(5,-5);    # nt 5 to fifth nt from the end
  $fq->substr_seq(-5)		# last 5 nts
  $fq->substr_seq(-100,50)  # 50 nts, starting 100 from end
  $fq->substr_seq(10,5,"AAAA") # replace 5 nts starting at pos 10 
    with 4 "A"s and replace corresponding qual values.
  ($fq1, $fq2) = $fq->substr_seq([5,100], [200,500]);
    # nts 5 to 100 and 200 to 500


=cut

sub substr_seq{
	my $self = shift;
	
	if(! ref $_[0] || @_ == 1){
		
		my $fa = $self->new; # clone

		my ($o, $l, $r) = ref $_[0] ? @{$_[0]} : @_;
					
		die __PACKAGE__."::substr_seq: Not enougth arguments\n" 
			unless defined ($o);
		
		# replace
		if(defined $r){
			$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
			substr($fa->{seq}, $o, $l, $r);
		}elsif(defined $l){
			$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
			$fa->seq( substr($fa->{seq}, $o, $l) );
		}else{
			$fa->desc_append(sprintf("SUBSTR:%d", $o));
			$fa->seq( substr($fa->{seq}, $o) );
		}
		return $fa;
	}else{
		my @new_fas;
		my $clone_c = 0;
		foreach (@_){
			$clone_c++;
			my $fa = $self->new; # clone
			$fa->id($fa->id.".$clone_c");
	
			my ($o, $l, $r) = @$_;
			
			die __PACKAGE__."::substr_seq: Not enougth arguments\n" 
				unless defined ($o);
			
			# replace
			if(defined $r){
				$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
				substr($fa->{seq}, $o, $l, $r);
			}elsif(defined $l){
				$fa->desc_append(sprintf("SUBSTR:%d,%d", $o, $l));
				$fa->seq( substr($fa->{seq}, $o, $l) );
			}else{
				$fa->desc_append(sprintf("SUBSTR:%d", $o));
				$fa->seq( substr($fa->{seq}, $o) );
			}
			push @new_fas, $fa; 
		}
		return @new_fas;
	}
	
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
	return $self->{desc};
}

=head2 string([LINEWIDTH])

Get entire sequence as FASTA string. Provide optional line width.

=cut

sub string{
	my ($self, $lw) = @_;
	my $s = $self->{seq};
	if($lw){
		$lw++;
		my $o=-1;
		substr($s, $o, 0, "\n") while (($o+=$lw) < length($s));
	}
	return sprintf("%s\n%s\n", $self->{seq_head}, $s);
}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=head2 seq_head

Get/Set the seq_head. Also updates id and desc.

=cut


sub seq_head{
	my ($self, $seq_head) = @_;
	if ($seq_head){
		$self->{seq_head} = $seq_head;
		# make sure, head has leading '>'
		substr($self->{seq_head},0,0,'>') unless $self->{seq_head} =~ /^>/;
		# reset id cache if head is changed
		my($id,$desc) = shift =~ m/
			(?:>?(\S+))			# id, >? for records
			(?:\s([^\n]+))?		# desc, optional
		/xs;	
		$self->{id} = $id; 
		$self->{desc} = $desc;
	};
	return $self->{seq_head};
}

=head2 seq

Get/Set the seq.

=cut

sub seq{
	my ($self, $seq) = @_;
	if($seq){
		$seq =~tr/\n//d;
		$self->{seq} = $seq 
	};
	return $self->{seq};
}

=head2 byte_offset

Get/Set the byte_offset.

=cut

sub byte_offset{
	my ($self, $byte_offset) = @_;
	$self->{byte_offset} = $byte_offset if defined($byte_offset);
	return $self->{byte_offset};
}


=head2 id

Get/set the seqs id. Also updates seq_head.

=cut

sub id{
	my ($self, $id) = @_;
	if(defined $id){
		$self->{id} = $id;
		# reset seq_head
		# make sure, desc is cached
		$self->{seq_head} = $self->{desc}
			? '>'.$self->{id}.' '.$self->{desc}
			: '>'.$self->{id};
	}
	return $self->{id};
}

=head2 desc

Get/set the seqs description, also updates seq_head.

=cut

sub desc{
	my ($self, $desc) = @_;
	if(defined $desc){
		$self->{desc} = $desc;
		# reset seq_head
		# make sure, desc is cached
		$self->{seq_head} = $self->{desc}
			? '>'.$self->{id}.' '.$self->{desc}
			: '>'.$self->{id}; # desc is '' ->remove desc
	}
	return $self->{desc};
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



