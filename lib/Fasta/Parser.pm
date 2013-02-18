package Fasta::Parser;

# $Id$

use warnings;
use strict;

# preference libs in same folder over @INC
use lib '../';

use Fasta::Seq 0.02;

use List::Util;

our $VERSION = '0.08';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;



##------------------------------------------------------------------------##

=head1 NAME 

Fasta::Parser.pm

=head1 DESCRIPTION

Parser module for Fasta format files.

=head1 SYNOPSIS

TODO

=cut

=head1 CHANGELOG

=head2 0.08

=over

=item [BugFix] C<< $fp->sample_seqs >> backups and restores buffer.

=item [Feature] C<< $fp->guess_seq_length >>

=item [Feature] C<< $fp->guess_seq_count >>

=item [Feature] C<< $fp->sample_seqs >> returns randomly drawn seq objects 
 from parsed file.

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

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

=head2 sample_seqs

Sample reads from file. If used on pipe returns undef. Takes 
 one argument, the number of reads to sample, default 1000. If the file is
 smaller than 10MB, the file is read entirely, else reads from random 
 positions are read. Returns LIST of Fastq::Seq objects.
 
=cut

sub sample_seqs{
	my ($self, $n) = (@_, 1000);
	my $fh = $self->fh;
	my $i;
	my @seqs;
	
	return if $self->check_fh_is_pipe;

	my $buffer = $self->{_buffer}; # backup buffer state to restore after sampling

	my $file_pos = tell($fh);

	# can seek on file: sample random reads
	my $size = -s $fh;
	if($size < 10_000_000){
		$self->seek(0);
		my $fa;
		my @fas;
		push @fas, $fa while $fa = $self->next_seq();
		my @shuffled_fas = List::Util::shuffle(@fas);
		@seqs = @shuffled_fas > $n
			? @shuffled_fas[0..$n] 
			: @shuffled_fas; 
		use Data::Dumper;
		print Dumper(scalar @seqs);
	}else{
		$size -= $size/100; # reduce size by 1% to prevent sampling eof
	
		for($i=$n;$i<1;$i--){
			$self->seek(int(rand($size))); # jump to random pos
			local $/ = "\n>";
			scalar <$fh>; # make sure to get a record start
			my $fa = $self->next_seq(); # get next record
			
			if($fa){
				push @seqs, $fa;
			}else{	# eof
				$i++; # do one more iteration
			}; 
		}
	}

	# restore file handle
	$self->seek($file_pos);
	
	$self->{_buffer} = $buffer; # restore buffer
	
	return @seqs;
}


=head2 guess_seq_length

Reads N sequences, randomly sampled if input is file and returns rounded 
 READLENGTH,STDDEV or undef on pipe, failure or empty file. Provide N as 
 the first parameter, default 1000.

=cut

sub guess_seq_length{
	my ($self, $n) = (@_, 1000);
	my $fh = $self->fh;
	my $i;
	
	my @sample_seqs = $self->sample_seqs($n);
	
	return undef unless @sample_seqs;
	
	my $l_total;
	my @lengths = map{my $l = length($_->seq); $l_total+=$l; $l}@sample_seqs;
	
	# empty file
	return undef unless $l_total;
	
	my $mean_l = $l_total/@sample_seqs;
	my $stddev = _stddev(\@lengths, $mean_l);
	# round mean
	return (int($mean_l + 0.5), int($stddev + 0.5)); 
}


=head2 guess_seq_count

Reads up to n reads from the current position of the file and estimates the
 mean size in bytes per FASTA record. Extrapolates the number sequences to 
 match the file size and returns the thus approximated total number of 
 sequences. Returns undef on STDIN.

=cut

sub guess_seq_count{
	my ($self, $n) = (@_, 1000);
	my $fh = $self->fh;
	# pipe
	return undef if $self->check_fh_is_pipe;

	my $size;
	my $file_size = -s $fh;
	# empty file
	return 0 unless $file_size;
	my @sample_seqs = $self->sample_seqs($n);
	$size+= length($_->string) for @sample_seqs;
	#print $median_size;
	return int(($file_size/$size* @sample_seqs)+0.5);
}

=head2 check_fh_is_pipe

Returns TRUE (1) if filehandle is associated with a pipe, FALSE ('') if it is 
 associated with a real file.

=cut

sub check_fh_is_pipe{
	my $self = shift;
	return $self->{_is_pipe};
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



