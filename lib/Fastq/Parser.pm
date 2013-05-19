package Fastq::Parser;

# $Id$

use warnings;
use strict;

use IO::File;
use IO::Uncompress::Gunzip;
use List::Util;

# preference libs in same folder over @INC
use lib '../';

use Fastq::Seq 0.10;


our $VERSION = '0.09';
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

=head2 0.09

=over

=item [Change] C<< $fp->next_seq(find_record => 1) >> replaces C<< $fp->next_seq(1) >>.
 Adjusted C<< $fp->sample_seq() >>.

=item [Change] Removed C<< $fp->next_raw_seq >>. Use C<< $fp->next_seq->string >>
 instead.

=back

=head2 0.08

=over

=item [BugFix] Changed C<< $fh->tell >> to C<< tell($fh) >>. 
 C<< IO::Handle->tell >> is only supported in latest module version.

=item [Change] STDIN is not dupped anymore. Dupped STDIN prevents subsequent
 reading of STDIN in main.

=item [Feature] << $fp->check_format >> now reads and B<unreads> the first
 char form input to determine format. Unreading makes it safe to use on 
 STDIN without using up stuff from the stream.

=item [Change] Removed << $fp->check_fh_is_pipe >>. Use << $fp->is_fh('PIPE') >> 
instead.
 
=item [Feature] << $fp->is_fh() >> allows to tell the type of input stream.

=item [Refactoring] Type of input stream is determined more sophistically. 

=item [Change] Uses <&STDIN dup instead of STDIN by default.

=back

=cut

=head2 0.07

=over

=item [BugFix] Filehandle to pipe can be -p or -t.

=item [Rename] C<< $fp->sample_seq >> to C<< $fp->sample_seqs >> 

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

=item empty file currently throws "NOT FASTA" exception

=back

=cut

##------------------------------------------------------------------------##

=head1 Class METHODS

=cut


##------------------------------------------------------------------------##

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
		fh => undef,
		file => undef,
		phred_offset => undef,
		mode => '<',
		_buffer => [],
		@_	# overwrite defaults
	};

	bless $self, $class;

	my $fh;
	# open file in read/write mode
	if($self->{file}){
		if(-T $self->{file}){
			open ($fh, $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
		}else{ # assume gzip input
			$fh = IO::Uncompress::Gunzip->new($self->{file});
		}
	}

	if($fh){
		$self->fh($fh);
	}else{
		#open(my $fh, "<&STDIN") or die $!;
		$self->fh(\*STDIN); # cannot dup, 'cause than ungetc does not work anymore!
	}

	return $self;
}

sub DESTROY{
	# just to be sure :D
	my $self = shift;
	close $self->fh unless $self->is_fh('PIPE');
}

=head1 Object METHODS

=cut

=head2 next_seq()

Loop through fastq file and return next 'Fastq::Seq' object. Without 
 parameter, assumes consistent, "four line wise" run through file,
 therefore it will fail if the filehandle is manually set to arbitrary 
 position. Returns undef on eof.

To perform a format checks on the seq object see the Fastq::Seq documention,
 C<< Fastq::Seq->CheckFormat() >>.

use << fp->next_seq('find_record => 1') >> to find the actual start of the 
 next valid record before retrieving it. This behaviour is useful in 
 combination with any previous C<< $fp->seek >> actions.


  $first_seq = $fp->next_seq();

  $fp->seek(-1000,2); # go 1000 bytes back from eof
  $close_to_last_seq = $fp->next_seq(find_record => 1);
  
  Fastq::Seq->CheckFormat(1); # global setting
  $checked_seq = $fp->next_seq();

=cut

sub next_seq{
	my $self = shift;
	my %p = (
		check_format => 0, 
		@_
	);
	
	
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
		
		if(exists $p{find_record} && $p{find_record}){
			# safe record start from any point in file
			my $lines = 0;
			until($nh =~ /^@/ && $qh =~ /^\+/){
				($nh,$ns,$qh,$qs) = ($ns,$qh,$qs, scalar <$fh>);
				
				# eof
				unless (defined $qs){
					# dont know if this makes sense
					# seek($fh,0,0); # reset to file start
					return;
				};
	
				# corrupt file
				die	sprintf("%s: %s, line %d, %s",(caller 0)[3],$self->{file}, $fh->input_line_number, 
					"Couldn't find record start within last 5 lines, possibly corrupted file or wrong format")	
					if $lines > 4;
				
				$lines++;
			}
		}
		
		# return fastq seq object
		return Fastq::Seq->new(
			$nh,$ns,$qh,$qs,
			phred_offset => $self->{phred_offset},
		);
	}
	
	return;
}


=head2 check_format

Takes a peek at the first entry in the file and checks wether the format of 
 the input looks like FASTQ (leading @). Returns the Parser object on 
 success, undef on failure. Does not modify the input stream, therefore
 can be used on STDIN safely.

NOTE: It only works at the start of the input. This means for pipes, use it
 before you start reading, on files use it either in the beginning or seek
 to the start to perform the check.

=cut

sub check_format{
	my ($self) = @_;
	my $fh = $self->fh;
	die sprintf("%s: %s",(caller 0)[3],"Format checking only works at the start of the file") 
		if tell($fh);
	my $c =$fh->getc(); # read first char
	# unread first char
	$self->is_fh('GZIP') 
		? $fh->ungetc($c)		# IO::Uncompress::Gunzip->ungetc pushes back string 
		: $fh->ungetc(ord($c)); # IO::File->ungetc pushes back char by ordinal

	return $c eq '@' ? $self : undef;
}


# DEPRECATED
#=head2 next_raw_seq
#
#Like C<< $fp->next_seq >> but returns next seq as raw string instead of 
# Fastq::Seq object. Setting the first parameter to TRUE also makes it 
# record safe.
#
#=cut

sub next_raw_seq{
	die "Use of next_raw_seq is deprecated, use something like next->seq->string instead";
#	my ($self, $safe) = (@_, 0);
#	my $fh = $self->{fh};
#	
#	if(@{$self->{_buffer}}){
#		# return fastq seq string
#		return join("",	splice(@{$self->{_buffer}}, 0, 4));
#	}
#	
#	if(
#		defined(my $nh = <$fh>) &&
#		defined(my $ns = <$fh>) &&
#		defined(my $qh = <$fh>) &&
#		defined(my $qs = <$fh>)
#	){
#		
#		if($safe){
#			# safe record start
#			my $lines = 0;
#			until($nh =~ /^@/ && $qh =~ /^\+/){
#				$nh = $ns;
#				$ns = $qh;
#				$qh = $qs;
#				$qs = <$fh>;
#				
#				# eof
#				unless (defined $qs){
#					seek($fh,0,0); # reset to file start
#					return;
#				};
#								# corrupt file
#				die	sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, "Couldn't find record start within next 5 lines, possibly corrupted file or wrong format")
#					if $lines > 4;
#				
#				$lines++;
#			}
#		}
#		
#		# return fastq seq string
#		return $nh.$ns.$qh.$qs;
#	}
#	
#	return;
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

=head2 sample_seqs

Sample reads from file. If used on pipe, returns N reads from the current
 position, while keeping them in the buffer for further processing. Takes 
 one argument, the number of reads to sample, default 1000. If the file is
 smaller than 10MB, the file is read entirely, else reads from random 
 positions are read. Returns LIST of Fastq::Seq objects.
 
NOTE: The set of samples reads is entirely kept in memory, therefore this
 method is not suited to sample large amount of reads from a large file.

=cut

sub sample_seqs{
	my ($self, $n) = (@_, 100);
	my $fh = $self->fh;
	my $i;
	my @reads;
	
	if($self->is_fh('PIPE') or $self->is_fh('GZIP')){
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
		my $size = $self->file_size;

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
				my $fq = $self->next_seq(find_record => 1); # get next reads with safe start
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
	
	my @sample_seq = $self->sample_seqs($n);
	
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
	
	my @sample_seq = $self->sample_seqs($n);
	
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
	# pipe
	return undef if $self->is_fh('PIPE');

	my $size;
	my $file_size = $self->file_size();
	# empty file
	return 0 unless $file_size;
	
	my @sample_seqs = $self->sample_seqs($n);
	return undef unless @sample_seqs;
	
	$size+= length($_->string) for @sample_seqs;
	#print $median_size;
	return int(($file_size/$size* @sample_seqs)+0.5);
}

=head2 is_fh

Determine the type of the filehandle. Without parameter, returns 0 for 
 handle to FILE, 1 for PIPE, and 2 for a handle to a SCALAR.
Alternatively you can provide the name of the type or the corresponding
 INT as single parameter. In these cases, the methods returns 1, if the
 type is matched, 0 otherwise.

  $fp->is_fh();  # 0,1 or 2
  $fp->is_fh('FILE') # 1 if filehandle is a handle to a FILE, 0 otherwise
  $fp->is_fh(0) # 1 if filehandle is a handle to a FILE, 0 otherwise

=cut

sub is_fh{
	my ($self, $type) = @_;
	
	my %type = (
		'FILE' => 0,
		'PIPE' => 1,
		'SCALAR' => 2,
		'GZIP' => 3,
		0 => 0,
		1 => 1,
		2 => 2,
		3 => 3,
	);
	
	if(defined $type){
		die sprintf("%s: %s",(caller 0)[3],"unknown type $type") 
			unless exists $type{$type};
		
		return $type{$type} == $self->{_is_fh} ? 1 : 0;
	}

	return $self->{_is_fh};
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
		if((ref $fh) =~ m/Gunzip/){
			$self->{_is_fh} = 3;
		}elsif(-f $fh){
			$self->{_is_fh} = 0;
		}elsif(-p $fh or  -t $fh){
			$self->{_is_fh} = 1;
		}elsif(ref $fh eq 'SCALAR'){
			$self->{_is_fh} = 2;
		}else{
			die sprintf("%s: %s",(caller 0)[3],"Unknown filehandle type");
		}
		$self->{fh} = $fh;
	}
	
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


=head2 file_size

Get the (uncompressed) size in kb of the underling data source.
Returns undef on pipes.

=cut

sub file_size{
	my ($self) = @_;
	if($self->is_fh('GZIP')){
		my @raw = `gzip --list $self->{file}`;
		my $size = ( split " ", $raw[1] )[1];
		return $size;
	}else{
		return -s $self->{fh};
	}
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



