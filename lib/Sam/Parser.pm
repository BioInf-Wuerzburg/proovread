package Sam::Parser;

use warnings;
use strict;

# $Id: Parser.pm 605 2012-12-09 12:03:56Z dumps $

# preference libs in same folder over @INC
use lib '../';

use Sam::Alignment qw(:flags);


our $VERSION = '0.11';




=head1 NAME 

Sam::Parser.pm

=head1 DESCRIPTION

Parser module for SAM format files.

=cut

=head1 SYNOPSIS

  use Sam::Parser;
  use Sam::Alignment ':flags';
  
  my $sp = Sam::Parser->new( 
      # parse from a file
    file => '/some/sam/file.sam',
      #or from a file handle
    fh => \*SAM,
      # or even from a STRINGREF
    file => \$sam_string,
  );
  
  # print names of all reference sequences from header
  while(%h = $sp->next_header_line('@SQ')){
  	print $h{'SN'}."\n";
  }
  
  # parser for file handle and with customized is routine 
  # read starts with 'C'
  my $sp = Sam::Parser->new( 
    fh => \*SAM,
    is => sub{ substr($_[0]->seq, 0, 1) eq 'C' }
  );	
  
  # print read ids of all reads with bad quality
  while( my $aln = $sp->next_aln() ){
    print $aln->qname() if $aln->is_bad_quality();
  }
  
  # seek the begin of the alignments for reparsing
  $sp->seek_alignment_section();
  
  # reset the 'is' routine
  $sp->is(MAPPED_BOTH);
  
  # print sequences of read pairs with both reads mapped
  while( my ($aln1, $aln2) = $sp->next_pair() ){
    print $aln1->seq().", ".$aln2->seq()."\n";
  }

=cut

=head1 CHANGELOG

=over

=head2 0.11

=over

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=item 0.10 [Thomas Hackl 2012-10-25]

Renamed Method. C<< $sp->append_tell() >> is now simply C<< $sp->tell() >>.

=item 0.09 [Thomas Hackl 2012-10-25]

Added C<< $sp->seek >> method.

=item 0.08 [Thomas Hackl 2012-10-25]

C<< $sp->next_header_line >> in LIST context now stores the raw line under 
 the key C<raw>.

=item 0.07 [Thomas Hackl 2012-10-24]

Performance upgrade. C<< _is >> tests are now not always performed, but only
 after testing if routine is defined.

=item 0.06 [Thomas Hackl 2012-10-24]

Added context awareness to C<< $sp->next_header_line >>. Returns now either
 the entire line as string or a LIST of splited tag => value pairs.

=item 0.05

Added C<< $sp->aln_by_pos >> to retrieve aln based on file position.

=item 0.04

Added explicit read mode C<< < >> to the file open statement. This
 allows to use the parser on STRINGREFS as if they where files.

=item 0.03

Renamed C<get_next_> to C<next_>, C<goto_> to C<seek_>.

Bugfixed eval part of C<is()> method. Object before method call was missing.

=item 0.02

Added C<_line_buffer> to object. This allows parsing without C<seek> and 
 therefore allows - besides reading from files - reading from STDIN.

Added header parsing support.

=item 0.01

Initial Parser module. Provides Constructor, generic accessor, individual
 and paired mode alignment parsing methods and conditional C<is> filtering 
 of alignments
 

=back

=cut

=head1 TODO

=over

=item Tests

=back

=cut

=head1 Constructor METHOD

=head2 new

Initialize a sam parser object. Takes parameters in key => value format. 

  fh => \*STDIN,
  file => undef,
  is => undef,
  mode => '<',   # read, 
                 # '+>': read+write (clobber file first)
                 # '+<': read+write (append)
                 # '>' : write (clobber file first)
                 # '>>': write (append)
=back

=cut

sub new{
	my $class = shift;
	
	my $self = {
		# defaults
		fh => \*STDIN,
		file => undef,
		is => undef,
		mode => '<',
		# overwrite defaults
		@_,
		# protected
		_line_buffer => undef,
		_aln_section => undef,
		_is => undef,
	};

	# open file in read/write mode
	if ($self->{file}){
		my $fh;
		open ( $fh , $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
		$self->{fh} = $fh;
	}

	bless $self, $class;
	
	# prepare _is test routine
	$self->is($self->{is}) if $self->{is};
	
	return $self;
	
}

sub DESTROY{
	# just to be sure :D
	my $self = shift;
	close $self->fh if $self->fh;
}








############################################################################


=head1 Object METHODS

=cut

=head2 next_aln

Loop through sam file and return next 'Sam::Alignment' object (meeting the 
 'is' criteria if specified).

=cut

sub next_aln{
	my %aln;
	my ($self, $seek) = @_;
	my $fh = $self->{fh};
	
	# process line buffer before looping
	if(my $sam = $self->{_line_buffer}){
		# clear buffer
		$self->{_line_buffer} = undef;
		
		$self->next_aln() if $sam =~ /^@/;
		if(!$self->{_aln_section}){
			$self->{_aln_section} = tell($fh)-length($sam);
		}
		# return sam aln object
		my $aln = Sam::Alignment->new($sam);
		return $aln if !$self->{_is} or &{$self->{_is}}($aln);
	}
	
	# loop if line buffer was empty or did not return
	while(<$fh>){
		# skip header
		if(!$self->{_aln_section}){
			next if /^@/;
			$self->{_aln_section} = tell($fh)-length($_);
		}

		# return sam aln object
		my $aln = Sam::Alignment->new($_);
		return $aln if !$self->{_is} or &{$self->{_is}}($aln);
	}
	#eof
	return;
}

=head2 next_pair / next_pair_single

Loop through sam file and return next 'Sam::Alignment' objects of a read 
 pair(matching the 'is' criteria if specified).

Which of the two methods you need depends on your data. The difference is 
 the handling of read pairs with only one mapped read ('half-mapped') by 
 different programs. If they are reported as single line you need to use 
 C<next_pair_single()> (e.g. segemehl with particular settings) else if
 they are reported in two lines use C<next_pair()>. Using the wrong 
 method will produce false results!

=over 12

=item next_pair

Use if half mapped pairs are reported as two lines/alignments in your sam. 
 Returns an list of the two alignment objects of a read pair.

=item next_pair_single

Use if half mapped pairs are reported as a single line/alignment in your sam. 
 Returns a list of the two alignment if both reads mapped and one alignment
 object if only one mapped.

=back

NOTE: C<next_pair(), next_pair_single()> can produce inconsistent 
 results if not used from the start of the file and if mixed with single 
 line operations like C<next_aln()> on the same parser instance.

=cut

sub next_pair{
	my ($self) = @_;
	my $fh = $self->{fh};
	
	# process line buffer before looping
	if(my $sam1 = $self->{_line_buffer}){
		# clear buffer
		$self->{_line_buffer} = undef;
		
		$self->next_pair() if $sam1 =~ /^@/;
		if(!$self->{_aln_section}){
			$self->{_aln_section} = tell($fh)-length($sam1);
		}
		# aln1
		my $aln1 = Sam::Alignment->new( $sam1 );
		# aln2
		my $aln2 = Sam::Alignment->new( scalar <$fh> );
		# return if property conditions are matched
		return ($aln1, $aln2) if !$self->{_is} || (&{$self->{_is}}($aln1) && &{$self->{_is}}($aln2));
	}

	# loop if line buffer was empty or did not return
	# get first aln of pair
	while(	my $sam1 = <$fh> ){
		# skip header
		if(!$self->{_aln_section}){
			next if $sam1 =~ /^@/;
			$self->{_aln_section} = tell($fh)-length($sam1);
		}
		
		# aln1
		my $aln1 = Sam::Alignment->new( $sam1 );
		# aln2
		my $aln2 = Sam::Alignment->new( scalar <$fh> );
		# return if property conditions are matched
		return ($aln1, $aln2) if !$self->{_is} || (&{$self->{_is}}($aln1) && &{$self->{_is}}($aln2));
	}
	#eof
	return;
}

sub next_pair_single{
	my ($self) = @_;
	my $fh = $self->{fh};
	
	# process line buffer before looping
	if(my $sam1 = $self->{_line_buffer}){
		# clear buffer
		$self->{_line_buffer} = undef;
		
		$self->next_pair_single() if $sam1 =~ /^@/;
		if(!$self->{_aln_section}){
			$self->{_aln_section} = tell($fh)-length($sam1);
		}
		
		# aln1
		my $aln1 = Sam::Alignment->new( $sam1 );			
		# aln2 ?
		unless($aln1->is(MAPPED_BOTH)){ # aln1 only
			return ($aln1, undef) if !$self->{_is} || &{$self->{_is}}($aln1)
		}else{ # aln2
			my $aln2 = Sam::Alignment->new( scalar <$fh> );
			return ($aln1, $aln2) if !$self->{_is} || (&{$self->{_is}}($aln1) && &{$self->{_is}}($aln2));
		}	
	}

	# loop if line buffer was empty or did not return
	# get first aln of pair
	while(	my $sam1 = <$fh> ){
		# skip header
		if(!$self->{_aln_section}){
			next if $sam1 =~ /^@/;
			$self->{_aln_section} = tell($fh)-length($sam1);
		}
		
		# aln1
		my $aln1 = Sam::Alignment->new( $sam1 );			
		# aln2 ?
		unless($aln1->is(MAPPED_BOTH)){ # aln1 only
			return ($aln1, undef) if !$self->{_is} ||  &{$self->{_is}}($aln1)
		}else{ # aln2
			my $aln2 = Sam::Alignment->new( scalar <$fh> );
			return ($aln1, $aln2) if !$self->{_is} ||  (&{$self->{_is}}($aln1) && &{$self->{_is}}($aln2));
		}
	}
	#eof
	return;	
}


=head2 aln_by_pos

Get an alignment, based on a byte offset position, return 'Sam::Alignment' 
 object (if meeting the 'is' criteria if specified) or undef. Also resets 
 the filehandle for C<< next_... >> methods.

  $aln = $sp->aln_by_pos($pos);

=cut

sub aln_by_pos{
	my ($self, $seek) = @_;
	my $fh = $self->{fh};
	
	seek($fh, $seek, 0) || return undef;
	
	my $l = scalar <$fh>;
	return if $l =~ m/^@/; # header section

	# return sam aln object
	my $aln = Sam::Alignment->new($l);
	return $aln if !$self->{_is} || &{$self->{_is}}($aln);
}

=head2 next_header_line

Parse linewise through sam file header information. The method returns the 
 entire line if used in SCALAR context, in LIST context a tag => value list 
 corresponding to the TAG:VALUE format of sam header lines, for convenience 
 each @CO comment line is returned as CO => <comment> and can be written 
 directly to a hash just like the other lines. In LIST context, the C<raw> 
 key contains the entire line. 
 
To retrieve a specific header line, provide the corresponding header 
 subsection key. If no key is given, any next header line is returned. 
 Returns FALSE if no (more) matching header lines are in the file. 

The following 
 subsection keys can but don't need to be present in a standard sam file.

  #key   meaning (mandatory TAGs of subsection entry)
  @HD => header line (VN)
  @SQ => reference sequence directory (SN, LN)
  @RG => read groups (ID)
  @PG => program (ID)
  @CO => comments
  
  # print names of all reference sequences from header
  while(%h = $sp->next_header_line('@SQ')){
  	print $h{'SN'}."\n";
  }
  

=cut

sub next_header_line{
	my ($self, $search_tag) = (@_, '@');
	my $fh = $self->{fh};
	
	# process line buffer before looping
	if(my $sam = $self->{_line_buffer}){
		unless($sam =~ /^@/){ # not header section
			return; # and keep line buffer filled
		}else{ # clear buffer
			$self->{_line_buffer} = undef;
		}
		
		if (my ($tag, $content) = $sam =~ /^(\@?(?:$search_tag)\w{0,2})\s(.*)/){
			if(wantarray){
				return $tag eq '@CO' 
					? (CO => $content, raw => $sam, tag => $tag) 
					: (split(/[:\t]/, $content), raw => $sam, tag => $tag);
			}else{
				return $sam;
			}
		}
	}
	
	# loop if line buffer was empty or did not return
	# get next header line
	while(	my $sam = <$fh> ){
		unless($sam =~ /^@/){ # end of header section
			$self->{_line_buffer} = $sam;
			return; 
											#  /^\@?($search_tag\w{0,2})\s(.*)/
		}elsif (my ($tag, $content) = $sam =~ /^(\@?(?:$search_tag)\w{0,2})\s(.*)/){
			if(wantarray){
				return $tag eq '@CO' 
					? (CO => $content, raw => $sam, tag => $tag) 
					: (split(/[:\t]/, $content), raw => $sam, tag => $tag);
			}else{
				return $sam;
			}
		}
		next;
	}
	return;
}

=head2 seek_alignment_section

Reset the file handle to the start of the alignment section. Resets 
 C<next_aln(), next_pair()> to the first aln in the sam file.

NOTE: this operation does only work on real files, not on STDIN.

=cut

sub seek_alignment_section{
	my ($self) = @_;
	die "".((caller 0)[3]).": Filehandle is a pipe, operation requires real file!" unless -f $self->fh;
	$self->{_line_buffer} = undef;
	seek($self->fh, $self->{_aln_section} ? $self->{_aln_section} : 0 ,0);
	return $self;
}

=head2 seek_header_section

Reset the file handle to the start of the header section (beginning of the
 file). Allows you to reread the header information.

NOTE: this operation does only work on real files, not on STDIN.

=cut


sub seek_header_section{
	my ($self) = @_;
	die "".((caller 0)[3]).": Filehandle is a pipe, operation requires real file!" unless -f $self->fh;
	$self->{_line_buffer} = undef;
	seek($self->fh, 0,0);
	return $self;
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


=head2 append_aln

Append an alignment to the file, provided as object or string. Returns the
 byte offset position in the file.

NOTE: In case a string is provided, make sure it contains trailing newline 
 since no further test is performed.

=cut

sub append_aln{
	my ($self, $aln) = @_;
	my $pos = tell($self->{fh});
	print {$self->{fh}} ref $aln ? $aln->raw : $aln;
	return $pos;
}


=head2 tell

Return the byte offset of the current append filehandle position

=cut

sub tell{
	my ($self) = @_;
	return tell($self->{fh});
}

=head2 append_tell

DEPRECATED: use C<< $fp->tell() >>

Return the byte offset of the current append filehandle position

=cut

sub append_tell{
	shift->tell(@_)
}


############################################################################

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


=head2 is

Get/Set conditions that determine which alignments are returned by the 
 next methods. Takes either a reference to a list of property bitmasks 
 or a code reference to a customized test which returns 1 and 0 respectively.
 To explicitly deactivate testing, provide a value that evaluates to FALSE.
 For details on bitmasks see L<Sam::Alignment>.
 
The test routine is executed with the parameters C<$parser_obj, $aln_obj> 
 and for C<next_pair()> additionally with C< $aln_obj2 >.

  # parser returning only BAD_QUALITY alns
  my $sp = Sam::Parser->new(
  	is => [Sam::Alignment->BAD_QUALITY]
  );
  
  # customized parser that only returns reads with a GC content > 70%.
  my $sp = Sam::Parser->new(
  	is => sub{
  	my ($self, $aln) = @_;
  	return ($aln->seq =~ tr/GC//) / length($aln->seq) > .7 ? 1 : 0;
  })
  
  # deactivate testing
  my $sp->is(0);
  
  
  
=cut

sub is{
	my ($self, $is) = @_;
	if(@_== 2){
		unless($is){
			$self->{_is} = undef;
		}elsif(ref($is) eq 'ARRAY'){
			$self->{_is} = eval 'sub{$_[0]->is('.join(', ', @$is).')}';
		}elsif(ref($is) eq 'CODE'){
			$self->{_is} = $is;
		}else{
			die (((caller 0)[3])." neither ARRAY nor CODE reference given!\n");
		}	
	}
	return $self->{_is};
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;



