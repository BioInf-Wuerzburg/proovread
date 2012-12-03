package Fastq::Parser;

use warnings;
use strict;
use Fastq::Seq 0.03;

our $VERSION = '0.06';

=head1 NAME 

Fastq::Parser.pm

=head1 DESCRIPTION

Parser module for FASTQ format files.

=head1 SYNOPSIS

TODO

=cut

=head1 CHANGELOG

=over

=item 0.06 [Thomas Hackl 2012-10-29]

Bugfix. << $fp->seek >> now clears _buffer.

=item 0.05 [Thomas Hackl 2012-10-28]

Added FASTQ write support, including byte offset indexing.

=item 0.04 [Thomas Hackl 2012-10-01]

Added C<< $fp->check_format >>. Determines whether input looks like FASTQ.

=item 0.03 [Thomas Hackl 2012-09-02]

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

=head1 Constructor METHOD

=head2 new

Initialize a fastq parser object. Takes parameters in key => value format. 
 Parameters are:

  fh => \*STDIN,
  file => undef,
  phred_offset => 64.
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
		fh => \*STDIN,
		file => undef,
		phred_offset => 64,
		mode => '<',
		_buffer => [],
		@_	# overwrite defaults
	};

	# open file in read/write mode
	if ($self->{file}){
		my $fh;
		open ( $fh , $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
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

=head1 Object METHODS

=cut

=head2 next_seq

Loop through fastq file and return next 'Fastq::Seq' object.

=cut


sub next_seq{
	my ($self) = @_;
	my $fh = $self->{fh};
	
	while(
		defined(my $nh = @{$self->{_buffer}} ? shift @{$self->{_buffer}} : <$fh>) &&
		defined(my $ns = <$fh>) &&
		defined(my $qh = <$fh>) &&
		defined(my $qs = <$fh>)
	){
		# return fastq seq object
		return Fastq::Seq->new(
			$nh,$ns,$qh,$qs,
			phred_offset => $self->{phred_offset}
		);
	}
}




=head1 Object METHODS

=cut

=head2 check_format

Check wether the format of the input looks like FASTQ (leading @). Returns
 the Parser object on success, undef on failure.

=cut

sub check_format{
	my $self = shift;
	my $fh = $self->fh;
	# read a line from/to the buffer
	$self->{_buffer} = [scalar <$fh>] unless @{$self->{_buffer}};
	my $l =	$self->{_buffer}[0];
	return $l =~ /^\@/ ? $self : undef;
}


=head2 next_raw_seq

Loop through fastq file and return next raw fastq sequence string. 

=cut

sub next_raw_seq{
	my ($self) = @_;
	my $fh = $self->{fh};
	
	while(
		defined(my $nh = @{$self->{_buffer}} ? shift @{$self->{_buffer}} : <$fh>) &&
		defined(my $ns = <$fh>) &&
		defined(my $qh = <$fh>) &&
		defined(my $qs = <$fh>)
	){
		# return fastq seq object
		return $nh.$ns.$qh.$qs;
	}
}

=head2 seek

Set the filehandle to the specified byte offset. Takes two
optional arguments "POSITION" (0), "WHENCE" (0), see perl "seek" for more.
Returns 'true' on success.

NOTE: this operation does only work on real files, not on STDIN.

=cut

sub seek{
	my ($self, $offset, $whence) = (@_, 0, 0);
	$self->{_buffer} = [];
	return seek($self->fh, $offset, $whence);
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

=head2 phred_offset

Get/Set phred quality offset. Default 64.

=cut

sub phred_offset{
	my ($self, $phred_offset) = @_;
	$self->{phred_offset} = $phred_offset if defined($phred_offset);
	return $self->{phred_offset};
}




=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



