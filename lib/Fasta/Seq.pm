package Fasta::Seq;

use warnings;
use strict;

use Verbose;

our $VERSION = '0.04';


##------------------------------------------------------------------------##

=head1 NAME 

Fasta::Seq.pm

=head1 DESCRIPTION

Class for handling FASTA sequences.

=head1 SYNOPSIS

=cut


=head1 CHANGELOG

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
	my $class = shift;
	my $self;
	
	if(@_%2){ # input is string to split
		my %self;
		@self{'id','desc', 'seq'} = shift =~ m/
			(?:>?(\S+))			# id, >? for records
			(?:\s([^\n]+))?\n	# desc, optional
			(.+)				# seq
		/xs;					# . includes \n
		$self = {
			byte_offset => undef,
			%self,
			@_	# overwrite defaults
		};
		# TODO: seq_head -> id -> id -> seq_head...
		$self->{seq_head} = '>'.$self->{id};
		$self->{seq_head}.= ' '.$self->{desc} if $self->{desc};

	}else{
		$self = {
			@_	# overwrite defaults
		};
		
		# make sure, id has not leading '>'
		$self->{id} =~ s/^>// if $self->{id};
		# make sure there is a seq_head entry
		if(!$self->{seq_head} && $self->{id}){
			$self->{seq_head} = '>'.$self->{id};
			$self->{seq_head}.= ' '.$self->{desc} if $self->{desc};
		}elsif(!$self->{id} && $self->{seq_head}){
			my($id,$desc) = shift =~ m/
				(?:>?(\S+))			# id, >? for records
				(?:\s([^\n]+))?		# desc, optional
			/xs;	
			$self->{id} = $id; 
			$self->{desc} = $desc;
		}else{
			$V->exit("Either seq_head or id required");
		}
	
	}

	chomp($self->{seq});
	chomp($self->{seq_head});
	$self->{seq} =~ tr/\n//d; 	# remove all newlines from seq

	return bless $self, $class;
}







##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2 clone

Create a clone of a seq object. Comes in handy of you want to modify the 
 object but still keep the original.
 
NOTE: This is not deep cloning but should work for this Class.

=cut

sub clone{
	my $self = shift;
	return bless ({%$self}, ref $self);
}


=head2 base_content

=cut

sub base_content{
	my $self = shift;
	my $patt = shift;
	__PACKAGE__->Add_base_content_scan($patt) unless exists $Base_content_scans->{$patt};
	return &{$Base_content_scans->{$patt}}($self->seq)
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
	$self->{seq} = $seq if $seq;
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

=head2 id

Get/set the seqs id. Also updates seq_head.

=cut

sub id{
	my ($self, $id) = @_;
	if(defined $id){
		$self->{id} = $id;
		# reset seq_head
		# make sure, desc is cached
		$self->{desc}
			? $self->{seq_head} = '>'.$self->{id}.' '.$self->{desc}
			: $self->{seq_head} = '>'.$self->{id};
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
		$self->{seq_head} = '>'.$self->{id}.' '.$self->{desc}
	}
	return $self->{desc};
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



