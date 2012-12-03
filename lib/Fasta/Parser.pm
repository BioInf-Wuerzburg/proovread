package Fasta::Parser;

use warnings;
use strict;
use Fasta::Seq 0.02;

our $VERSION = '0.07';





##------------------------------------------------------------------------##

=head1 NAME 

Fasta::Parser.pm

=head1 DESCRIPTION

Parser module for Fasta format files.

=head1 SYNOPSIS

TODO

=cut

=head1 CHANGELOG

=over

=item 0.07 [Thomas Hackl 2012/11/13]

Bugfix. C<< $fp->next_seq >> now correctly runs tell on 
 C<< $fp->fh >> and not latest selected FILEHANDLE (STDIN).

=item 0.06 [Thomas Hackl 2012/10/25]

Added C<< $fp->seek >> to set the filehandle to a specific position.

=item 0.05 [Thomas Hackl 2012-10-24]

Check format only reads the first char of the input and returns
 success it it matches '>'. Only reading the '>' makes _buffer 
 obsolete since C<< Fastq::Seq->new() >> does not require a leading '>'.

Removed the _buffer feature.

C<< $fap->next_seq() >> now also stores the byte offset of the 
 sequence returned in the object.

=item 0.04 [Thomas Hackl 2012/10/01]

Added C<< $fp->check_format >>. Determines whether input looks like FASTQ.

=item 0.03 [Thomas Hackl 2012/09/02]

Modified to use Fastq:Seq 0.03+

=item 0.02 [Thomas Hackl]

Added C<next_raw_seq()>.

=item 0.01 [Thomas Hackl]

Initial Parser module. Provides Constructor and generic accessor
 methods

=back

=cut

=head1 TODO

=over

=item Synopsis

=item Tests

=item read reads chunkwise

=back

=cut




##------------------------------------------------------------------------##

=head1 Constructor METHOD

=head2 new

Initialize a fasta parser object. Takes parameters in key => value format. 
 Parameters are:

  fh => undef,
  file => undef,
  mode => '<',   # read, 
                 # '+>': read+write (clobber file first)
                 # '+<': read+write (append)
                 # '>' : write (clobber file first)
                 # '>>': write (append)

=cut

sub new{
	my $class = shift;
	
	# defaults=
	my $self = {
		fh => undef,
		file => undef,
		mode => '<',
		@_	# overwrite defaults
	};

	# open file in read/write mode
	if ($self->{file}){
		my $fh;
		open ( $fh , $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
		$self->{fh} = $fh;
	}
	
	# create a STDIN alias for save parsing and closing
	unless($self->{fh}){
		open(my $fh, "<&STDIN") or die $!;
		$self->{fh} = $fh;
	}
	
	bless $self, $class;
	return $self;
}

sub DESTROY{
	# just to be sure :D
	my $self = shift;
	close $self->fh;
}



##------------------------------------------------------------------------##

=head1 Object METHODS

=cut

=head2 next_seq

Loop through fasta file and return next 'Fasta::Seq' object.

=cut


sub next_seq{
	my ($self) = @_;
	
	my $fh = $self->{fh};
	local $/ = "\n>";
	
	# return fasta seq object
	my $byte_offset = tell($fh);
	my $fa = <$fh>;
	return unless defined $fa;
	chomp($fa);
	return Fasta::Seq->new($fa, byte_offset => $byte_offset);
}




=head2 check_format

Check wether the format of the input looks like FASTA (leading >). 
 Returns the Parser object on success, undef on failure.

=cut

sub check_format{
	my ($self) = @_;
	my $fh = $self->fh;
	my $c;
	# read first char
	read($fh,$c,1,0);
	return $c eq '>' ? $self : undef;
}



=head2 append_seq

Append an sequence to the file, provided as object or string. Returns the
 byte offset position in the file.

NOTE: In case a string is provided, make sure it contains trailing newline 
 since no further test is performed.

=cut

sub append_seq{
	my ($self, $seq) = @_;
	my $pos = tell($self->{fh});
	print {$self->{fh}} ref $seq ? $seq->string : $seq;
	return $pos;
}


=head2 append_tell

Return the byte offset of the current append filehandle position

=cut

sub append_tell{
	my ($self) = @_;
	return tell($self->{fh});
}

=head2 seek

Set the filehandle to the specified byte offset. Takes two
optional arguments "POSITION" (0), "WHENCE" (0), see perl "seek" for more.
Returns 'true' on success.

NOTE: this operation does only work on real files, not on STDIN.

=cut

sub seek{
	my ($self, $offset, $whence) = (@_, 0, 0);
	return seek($self->fh, $offset, $whence);
}

##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 fh

Get/Set the file handle.

=cut

sub fh{
	my ($self, $fh) = @_;
	$self->{fh} = $fh if $fh;
	return $self->{fh};
}




##------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



