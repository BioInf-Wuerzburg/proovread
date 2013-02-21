package Fastq::Parser;

# $Id$

use warnings;
use strict;

# preference libs in same folder over @INC
use lib '../';

use Fastq::Seq 0.03;

use List::Util;

our $VERSION = '0.07';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (.+)/;



##------------------------------------------------------------------------##

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

=item [BugFix] C<< $fp->sample_seq >> backups and restores buffer.

=item [BugFix] C<< $fp->sample_seq >> did not tell file positiion on correct
 file handle.

=item [Change] C<< $fp->guess_phred_offset >> automatically sets 
 Parsers C<< $fp->{phred_offset} >> attribute.

=item [Change] Preference libs in same folder over @INC

=item [Change] C<< $fp->guess_seq_length >> now returns estimated mean
 READLENGTH, STDDEV. In case a file is provided, reads are randomly sampled,
 not read from the current file position as done with pipes. By default 
 reads 1000 reads.

=item [Feature] C<< $fp->next_seq >> and C<< $fp->next_raw_seq >> now take
 an optional boolean true which makes the methods search safely for the next 
 record entry, regardless of any arbitrary position the file handle currently
 might have. Useful after seek stuff.

=item [Feature] C<< $fp->check_fh_is_pipe >> checks whether filehandle is 
 associated with pipe or file.

=item [BugFix] Splicing reads from buffer had error.

=item [Change] Added auto Id, Date and Revision on svn checkin

=item [Feature] C<< $fp->guess_phred_offset >>

=item [Feature] C<< $fp->guess_seq_length >>

=item [Feature] C<< $fp->guess_seq_count >>

=item [Feature] C<< $fp->sample_seqs >> returns randomly drawn seq objects 
 from parsed file.

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

=back

=cut

=head1 Constructor METHOD

=head2 new

Initialize a Fastq::Parser object. Takes parameters in key => value format. 
 Parameters are:

  fh => \*STDIN,
  file => undef,
  phred_offset => 64.
  mode => '<',  # read, 
        # '+>'    read+write (clobber file first)
        # '+<'    read+write (append)
        # '>'     write (clobber file first)
        # '>>'    write (append)

=cut

sub new{
	my $class = shift;
	
	# defaults=
	my $self = {
		fh => \*STDIN,
		file => undef,
		phred_offset => undef,
		mode => '<',
		_buffer => [],
		_is_pipe => undef,
		@_	# overwrite defaults
	};

	# open file in read/write mode
	if ($self->{file}){
		my $fh;
		open ( $fh , $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
		$self->{fh} = $fh;
		$self->{_is_pipe} = '';
	}else{
		 $self->{_is_pipe} = -p $self->{fh} ? 1 : ''; 
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

=head2 next_seq()

Loop through fastq file and return next 'Fastq::Seq' object. Without 
 parameter, assumes consistent, "four line wise" run through file,
 therefore it will fail if the filehandle is manually set to arbitrary 
 position. 

Set the first parameter to TRUE, to perform an additional check to find the 
 actual start of the next valid record before retrieving it. This behaviour
 is useful in combination with any previous C<< $fp->seek >> actions.

Returns undef on eof.

  $first_seq = $fp->next_seq();

  $fp->seek(-1000,2); # go 1000 bytes back from eof
  $close_to_last_seq = $fp->next_seq(1);


=cut

sub next_seq{
	my ($self, $safe) = (@_, 0);
	
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
		
		if($safe){
			# safe record start
			my $lines = 0;
			until($nh =~ /^@/ && $qh =~ /^\+/){
				$nh = $ns;
				$ns = $qh;
				$qh = $qs;
				$qs = <$fh>;
				# eof
				return unless defined ($qs);
	
				# corrupt file
				die	sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, "Couldn't find record start within next 5 lines, possibly corrupted file or wrong format")
					if $lines > 4;
				
				$lines++;
			}
		}
		
		# return fastq seq object
		return Fastq::Seq->new(
			$nh,$ns,$qh,$qs,
			phred_offset => $self->{phred_offset}
		);
	}
	
	return;
}


=head2 check_format

Takes a peek at the first entry in the file and checks wether the format of 
 the input looks like FASTQ (leading @, + at start of third line). Returns 
 the Parser object on success, undef on failure.
 
NOTE: Use this directly after creating the Parser object, if other methods 
 are run on the Fastq::Parser object prior, it is likely to 
 fail.

=cut

sub check_format{
	my $self = shift;
	my $fh = $self->fh;
	# read a line from/to the buffer
	$self->{_buffer} = [my $nh = scalar <$fh>, scalar <$fh>, my $qh = scalar <$fh>, scalar <$fh>] unless @{$self->{_buffer}};
	if($nh =~ /^@/ and $qh =~ /^\+/){
		return $self;
	}
	return undef;
}


=head2 next_raw_seq

Like C<< $fp->next_seq >> but returns next seq as raw string instead of 
 Fastq::Seq object. Setting the first parameter to TRUE also makes it 
 record safe.

=cut

sub next_raw_seq{
	my ($self, $safe) = (@_, 0);
	my $fh = $self->{fh};
	
	if(@{$self->{_buffer}}){
		# return fastq seq string
		return join("",	splice(@{$self->{_buffer}}, 0, 4));
	}
	
	if(
		defined(my $nh = <$fh>) &&
		defined(my $ns = <$fh>) &&
		defined(my $qh = <$fh>) &&
		defined(my $qs = <$fh>)
	){
		
		if($safe){
			# safe record start
			my $lines = 0;
			until($nh =~ /^@/ && $qh =~ /^\+/){
				$nh = $ns;
				$ns = $qh;
				$qh = $qs;
				$qs = <$fh>;
				# eof
				return unless defined ($qs);
	
				# corrupt file
				die	sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, "Couldn't find record start within next 5 lines, possibly corrupted file or wrong format")
					if $lines > 4;
				
				$lines++;
			}
		}
		
		# return fastq seq string
		return $nh.$ns.$qh.$qs;
	}
	
	return;
}


=head2 seek

Set the filehandle to the specified byte offset. Takes two
optional arguments POSITION (default 0), WHENCE (default 0), see perl "seek" for more.
Returns 'true' on success.

NOTE: this operation only works on real files, not on STDIN.

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

=head2 sample_seq

Sample reads from file. If used on pipe, returns N reads from the current
 position, while keeping them in the buffer for further processing. Takes 
 one argument, the number of reads to sample, default 1000. If the file is
 smaller than 10MB, the file is read entirely, else reads from random 
 positions are read. Returns LIST of Fastq::Seq objects.
 
NOTE: The set of samples reads is entirely kept in memory, therefore this
 method is not suited to sample large amount of reads from a large file.

=cut

sub sample_seq{
	my ($self, $n) = (@_, 1000);
	my $fh = $self->fh;
	my $i;
	my @reads;
	
	if($self->check_fh_is_pipe){
		#cant seek on pipe: need to buffer and can only read head
		
		my $l;
		# read lines from buffer
		for($i=0; $i<$n && 	$i < (scalar @{$self->{_buffer}}-3)/4; $i++){
			push @reads,  Fastq::Seq->new(
				@{$self->{_buffer}}[($i*4)..($i*4+3)],
				phred_offset => $self->{phred_offset}
			);
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
			# return fastq seq object
			push @reads, Fastq::Seq->new(
				$nh,$ns,$qh,$qs,
				phred_offset => $self->{phred_offset}
			);
			$i++;		
		}
	}else{
		my $buffer = $self->{_buffer}; # backup buffer state to restore after sampling
		my $file_pos = tell($fh);
		# can seek on file: sample random reads
		my $size = -s $fh;

		if($size < 10_000_000){
			$self->seek(0);
			my $fq;
			my @fqs;
			push @fqs, $fq while $fq = $self->next_seq();
			my @shuffled_fqs = List::Util::shuffle(@fqs);
			@reads = @shuffled_fqs > $n
				? @shuffled_fqs[0..$n] 
				: @shuffled_fqs; 
		}else{

			$size -= $size/100; # reduce size by 1% to prevent sampling eof
	
			for($i=$n;$i;$i--){
				$self->seek(int(rand($size))); # jump to random pos
				my $fq = $self->next_seq(1); # get next reads with safe start
				if($fq){
					push @reads, $fq;
				}else{	# eof
					$i++; # do one more iteration
				}; 
			}
		}
		
		# restore file handle
		$self->seek($file_pos);
		$self->{_buffer} = $buffer; # restore buffer
	}
	
	return @reads;
}


=head2 guess_seq_length

Reads up to N reads, randomly sampled if input is file, from the current 
 position if input is pipe, and returns rounded (READLENGTH,STDDEV), or 
 undef on failure or empty file. Provide N as the first parameter, 
 default 1000.

=cut

sub guess_seq_length{
	my ($self, $n) = (@_, 1000);
	my $fh = $self->fh;
	my $i;
	
	my @sample_seq = $self->sample_seq($n);
	
	return undef unless @sample_seq;
	
	my $l_total;
	my @lengths = map{my $l = length($_->seq); $l_total+=$l; $l}@sample_seq;
	
	# empty file
	return undef unless $l_total;
	
	my $mean_l = $l_total/@sample_seq;
	my $stddev = _stddev(\@lengths, $mean_l);
	# round mean
	return (int($mean_l + 0.5), int($stddev + 0.5)); 
}

=head2 guess_phred_offset

Samples up to N reads from current data and checks the boundaries of the 
 quality range. Returns 33 for range (33,33+42), 64 for range (64,64+42)
 or undef in any other case. Provide N as the first 
 parameter, default 1000.
 
 NOTE: There is an intersection for both ranges (64 to 75). If all sampled 
 values lie within this intersection range the method will return undef, 
 since the offset cannot be determined with certainty.

=cut

sub guess_phred_offset{
	my ($self, $n) = (@_, 1000);
	my $fh = $self->fh;
	my $i;
	
	my @sample_seq = $self->sample_seq($n);
	
	return undef unless @sample_seq;
	
	my @quals;
	foreach my $fq(@sample_seq){
		push @quals, split(//, $fq->qual);
	}
	
	@quals = sort @quals;
	
	my $min = $quals[0];
	my $max = $quals[-1];
	
	# intersection => undetermined
	if( ord($min) >= 64 && ord($max) <= 33+42){
		$self->{phred_offset} = undef;
		return undef 
	
	# 33
	}elsif( ord($min) >= 33 && ord($max) <= 33+42 ){
		$self->phred_offset(33);
		return 33 
	
	# 64
	}elsif(	ord($min) >= 64 && ord($max) <= 64+42 ){
		$self->phred_offset(64);
		return 64 
	
	# unknown
	}else{
		$self->{phred_offset} = undef;
		return undef;
	}
	
}

=head2 guess_seq_count

Reads up to n reads from the current position of the FASTQ and estimates the
 mean size in bytes per FASTQ entry. Extrapolates the read number to match 
 the file size and returns the thus approximated total number of reads. 
 Returns undef on STDIN.

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
	
	my @sample_seqs = $self->sample_seq($n);
	return undef unless @sample_seqs;
	
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



=head1 Accessor METHODS

=cut

=head2 fh

Get/Set the file handle.

=cut

sub fh{
	my ($self, $fh) = @_;
	if($fh){
		$self->{fh} = $fh;
		$self->{_is_pipe} = -p $fh ? 1 : undef;	
	};
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



