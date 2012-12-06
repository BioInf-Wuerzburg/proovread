package Fastq::Parser;

# $Id$

use warnings;
use strict;
use Fastq::Seq 0.03;

our $VERSION = '0.07';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our $MODIFIED = '$Date$' =~ /Date: (.+)/;

=head1 NAME 

Fastq::Parser.pm

=head1 DESCRIPTION

Parser module for FASTQ format data. Reads from files, STDIN or other 
 perl filehandles, can be used to write FASTQ data as well. 
 
Most methods work on any kind of input, some methods, like seek will fail 
 on streams and work only on real files. 
 
Main feature of the module is to sequentially run through the data and 
 retrieve Fastq::Seq objects, one at a time for further processing.

=head1 SYNOPSIS

TODO

=cut

=head1 CHANGELOG

=head2 0.07

=over

=item [BugFix] Splicing reads from buffer had error.

=item [Change] Added auto Id, Date and Revision on svn checkin

=item [Feature] C<< $fp->guess_phred_offset >>

=item [Feature] C<< $fp->guess_read_length >>

=item [Feature] C<< $fp->guess_total_reads >>

=item [Change] C<< $fp->check_format >> now stores a complete entry (4 lines)

=back

=head2 0.06

[Bugfix] C<< $fp->seek >> now clears _buffer.

=head2 0.05

Added FASTQ write support, including byte offset indexing.

=head2 0.04

Added C<< $fp->check_format >>. Determines whether input looks like FASTQ.

=head2 0.03

Modified to use Fastq:Seq 0.03+

=head2 0.02

Added C<next_raw_seq()>.

=head2 0.01

Initial Parser module. Provides Constructor and generic accessor
 methods


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
	
	use Data::Dumper;
	
	if(@{$self->{_buffer}}){
		# return fastq seq object
		return Fastq::Seq->new(
			splice(@{$self->{_buffer}}, 0, 4),
			phred_offset => $self->{phred_offset}
		);
	}
	
	my $fh = $self->{fh};
	if(
		defined(my $nh = <$fh>) &&
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
	
	return;
}



=head2 check_format

Takes a peek at the next entry in the file and checks wether the format of 
 the input looks like FASTQ (leading @). Returns the Parser object on success, 
 undef on failure.

=cut

sub check_format{
	my $self = shift;
	my $fh = $self->fh;
	# read a line from/to the buffer
	$self->{_buffer} = [scalar <$fh>, scalar <$fh>, scalar <$fh>, scalar <$fh>] unless @{$self->{_buffer}};
	my $l =	$self->{_buffer}[0];
	return $l =~ /^\@/ ? $self : undef;
}


=head2 next_raw_seq

Loop through fastq file and return next raw fastq sequence string. 

=cut

sub next_raw_seq{
	my ($self) = @_;
	my $fh = $self->{fh};
	
	if(@{$self->{_buffer}}){
		# return fastq seq string
		return join("",	splice(@{$self->{_buffer}}, 4));
	}
	
	if(
		defined(my $nh = <$fh>) &&
		defined(my $ns = <$fh>) &&
		defined(my $qh = <$fh>) &&
		defined(my $qs = <$fh>)
	){
		# return fastq seq string
		return $nh.$ns.$qh.$qs;
	}
	
	return;
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

=head2 guess_read_length

Reads up to n reads from the current position of the FASTQ and returns either 
 READLENGTH if all n reads had identical lengths or shortcuts if to reads
 differ in length, returns undef. Provide n as the first parameter, default 
 100.

=cut

sub guess_read_length{
	my ($self, $n) = (@_, 100);
	my $fh = $self->fh;
	my $i;
	my $l;
	# read a lines from buffer
	for($i=0; $i<$n && $i < (@{$self->{_buffer}}-3)/4; $i++){
		my $ns = $self->{_buffer}[($i*4) +1];
		$l = length($ns)-1 unless $l; # unchomped
		return 0 unless $l == length($ns)-1;
	}
	
	# read new lines and add to buffer
	while(
		$i < $n && 
		defined (my $nh = scalar <$fh>) &&
		defined (my $ns = scalar <$fh>) &&
		defined (my $qh = scalar <$fh>) &&
		defined (my $qs = scalar <$fh>)
	){
		push @{$self->{_buffer}}, $nh, $ns, $qh, $qs;
		$l = length($ns)-1 unless $l; # unchomped
		return 0 unless $l == length($ns)-1;
		$i++;		
	}
	return $l;
}

=head2 guess_phred_offset

Reads up to n reads from the current position of the FASTQ and scans the quality 
 sequence. If a numeric value of a quality character below 64 is found, 
 returns 33, if a value above 33+43 (76) is found, returns 64. If neither
 condition is met within n reads, returns undef. Provide n as the first 
 parameter, default 100.

=cut

sub guess_phred_offset{
	my ($self, $n) = (@_, 100);
	my $fh = $self->fh;
	my $i;
	my $min;
	my $max;
	# read a lines from buffer
	for($i=0; $i<$n && $i < (@{$self->{_buffer}}-3)/4; $i++){
		my $qs = $self->{_buffer}[($i*4) +3];
		foreach(split(//, $qs)){
			return 33 if ord($_) < 64;
			return 64 if ord($_) > 76;
		}
	}
	
	# read new lines and add to buffer
	while(
		$i < $n && 
		defined (my $nh = scalar <$fh>) &&
		defined (my $ns = scalar <$fh>) &&
		defined (my $qh = scalar <$fh>) &&
		defined (my $qs = scalar <$fh>)
	){
		push @{$self->{_buffer}}, $nh, $ns, $qh, $qs;
		foreach(split(//, $qs)){
			return 33 if ord($_) < 64;
			return 64 if ord($_) > 76;
		}
		$i++;		
	}
}

=head2 guess_total_reads

Reads up to n reads from the current position of the FASTQ and estimates the
 mean size in bytes per FASTQ entry. Extrapolates the read number match the 
 file size and returns the thus approximated total number of reads. Returns
 undef on STDIN.

=cut

sub guess_total_reads{
	my ($self, $n) = (@_, 100);
	my $fh = $self->fh;
	my $i;
	my $size;
	my $file_size = -s $fh;
	return unless $file_size; # 
	# read a lines from buffer
	for($i=0; $i<$n && $i < (@{$self->{_buffer}}-3)/4; $i++){
		my @fq = @{$self->{_buffer}}[($i*4)..($i*4)+3];
		$size+=length(join("", @fq));
	}
	
	# read new lines and add to buffer
	while(
		$i < $n && 
		defined (my $nh = scalar <$fh>) &&
		defined (my $ns = scalar <$fh>) &&
		defined (my $qh = scalar <$fh>) &&
		defined (my $qs = scalar <$fh>)
	){
		push @{$self->{_buffer}}, $nh, $ns, $qh, $qs;
		$size+=length(join("", $nh, $ns, $qh, $qs));
		$i++;
	}
	
	return int(($file_size/$size)*$i);
	
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



