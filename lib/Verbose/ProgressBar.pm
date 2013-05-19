package Verbose::ProgressBar;

use warnings;
use strict;

# $Id$

# preference libs in same folder over @INC
use lib '../';

use Verbose;

use IO::File;
use IO::Uncompress::Gunzip;

our $VERSION = '0.08';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;

$|++;



##------------------------------------------------------------------------##

=head1 NAME 

Verbose::ProgressBar - Progress Visualisation.

=cut

=head1 DESCRIPTION

Visualizes progress in an self-updating pv like ProgressBar, inlcuding ETA
 calculation.

=cut

=head1 SYNOPSIS


  # Counter based progress
  my $iter = 10000;
  my $pg = Verbose::ProgressBar->new(
    size => $iter
  )
  
  my $i;
  for($i=0;$i++;$i<$iter){
  	$i++;
  	$pg->update($i)
  }
  $pg->finish($i);
  
  # File size based progress
  open(my FILE, "my.file");
  my $pg = Verbose::ProgressBar->new(
    size => $fh,
    line_width => 120,
  )
  # parse large file linewise
  my $c=0;
  while(<FILE>){
    # invoke update every 10000th line
    $pg->update unless $c%10000;
  };
  $pg->finish;


=cut

=head1 CHANGELOG

=head2 0.08

=over

=item [Feature] Added C<level> and C<report_level> feature to modify behaviour
 according to given verbose level.

=item [BugFix] TTS and ETA now also show number of days, if required.

=item [BugFix] TTS is now always printed in case finish is called.

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=over

=item 0.07

[BugFix] Update with undefined/null value created error.

Bugfix. Initial bar was one space short.

=item 0.06 [Thomas Hackl 2012-10-10]

C<< $pb->finish >> resets the ProgressBar so that it can be used again, e.g.
 if the same file is read again.

=item 0.05

After C<< $pb->finish() >> call, ETA field in output is replaced by TTS 
(total time spent).

=item 0.04

Bugfix, C<< $self->{value} >> was not set if exteral value is provided.

=item 0.03

Cooldown for C<< $pb->update() >> of one second. Prevents overhead if update
 is called very often.

=item 0.02

Added ETA. Added C<< $pb->finish() >> for finishing the output. Added 
 "endless" monitoring for unknown size.

=item 0.01

Initial Verbose::ProgressBar module. Provides constructor and 
 update methods. Computes and print a ProgressBar.

=back

=head1 Constructor METHOD

=cut


sub new{
	my $class = shift;
	
	# overwrite defaults
	my $self = {
		fh => \*STDERR,
		size => undef,		# number, file or fh
		line_width => 80,	#
		level => 1,
		report_level => 1, 
		_fh =>undef, 		# a monitored file handle
		_bins => undef,
		_bin => -1,
		_bin_time => undef,
		_start_time => undef,
		_time => 0,
		_template => undef,
		@_%2 ? (size => shift, @_) : @_	# overwrite + odd number of params -> first is message
	};
	
	# calculate max size if filehandle
	if(defined ($self->{size}) && ref $self->{size}){
		$self->{_fh} = $self->{size};
		$self->{size} = -s $self->{_fh};
	}

	# calculate number of bins 
	$self->{_bins} = $self->{line_width} - 23; # 7 PBval, 3 PBbar, 13 PBeta
	$self->{_template} = "\r%6s [%-".($self->{_bins}+1)."s]";

	bless $self, $class;
}


=head1 Class METHODS

=head2 update

NOTE: The update function has a cooldown of 1 second to prevent 
 computational overhead, still it can cost you some seconds. 
 If this is a problem, performance can be further increased using
 an external counter and invoking update only every x iterations.

  # parse large file linewise
  my $c=0;
  while(<FILE>){
  	 # invoke update every 10000th line
  	 $pg->update unless $c%10000;
  };

=cut

sub update{
	my ($self, $value) = @_;
	
	# short circuit if less than 1 second
	(my $time = time)-$self->{_time} || return;
	$self->{_start_time} || ($self->{_start_time} = $time);
	
	$self->{level} > $self->{report_level} && return;
	
	if (!defined ($value) && defined($self->{_fh})){
		$self->{value} = $self->{_fh}->tell;
	}elsif(!$value){
		$self->{value} = $value = 0;
	}else{
		$self->{value} = $value;
	}
	
	# progress bar with known size
	
	my $eta;
	if($self->{size}){
		# calculate current bin
		my $bin = int($self->{value} /$self->{size} * $self->{_bins});

		# return unless new bin
		return if $bin <= $self->{_bin} && $self->{_time} > -1;
		$self->{_bin} = $bin;
		
		# initialize
		# 
		unless ($bin){
			print {$self->{fh}} "    0  [".( " "x ($self->{line_width} - 22))."]";
			return;
		}
		
		# (re)calculate bin time
		$self->{_bin_time} = ($time-$self->{_start_time})/$bin;
		if($self->{_time} < 0){ # reset by finish
			my $ela_time = $time-$self->{_start_time};		
			$eta = $ela_time >= 86400 
				? sprintf(" TTS %dd %02d:%02d:%02d", (gmtime($ela_time))[7,2,1,0])
				: sprintf(" TTS %02d:%02d:%02d", (gmtime($ela_time))[2,1,0])
		}elsif(my $eta_time = ($self->{_bins}-$bin) * $self->{_bin_time}){
			$eta = $eta_time >= 86400 
				? sprintf(" ETA %dd %02d:%02d:%02d", (gmtime($eta_time))[7,2,1,0])
				: sprintf(" ETA %02d:%02d:%02d", (gmtime($eta_time))[2,1,0])
		}else{
			$eta = '';
		}
		printf {$self->{fh}} 
			($self->{_template}.$eta,
				Verbose->Humanize($self->{value}),		# PGval
				'='x($self->{_bin}).'>',		# PGbar
			);
		
	}else{
		# fill/reset bins in endless loop
		if($self->{_bin} < $self->{_bins}){
			$self->{_bin}++
		}else{
			$self->{_bin} = 0;
		}
		
		if(defined $value){
			$self->{value} = $value;
		}else{
			# use simple incrementor if no values are provided
			$self->{value}++ 
		}
		# print progress line
		printf {$self->{fh}} 
		($self->{_template},
			Verbose->Humanize($self->{value}),		# PGval
			' 'x($self->{_bin}-2).'<=>',		# PGbar
		);
	}
	
	# save time	
	$self->{_time} = $time;
}

=head2 finish

Call after monitoring is finished. Requires the same input as 
 C<update>, which means none in case a file handle is monitored or
 the final counter value in case an external counter is used.
Finishes the output.

Returns the final size and the total time required. Resets the bar, so it 
 can be used again.

=cut

sub finish{
	my ($self, $value) = @_;
	$self->{level} > $self->{report_level} && return;
	# run the last bin, if it hasn't been run. Can happen if resolution is
	# small und update on the final bin is not triggered
	$self->{_time} = -1; # disable short circuit
	$self->update($value);
	print {$self->{fh}} "\n";
	my @re = ($self->{value}, time()-$self->{_start_time});
	
	# reset
	$self->{_bin} = -1,
	$self->{_bin_time} = undef;
	$self->{_start_time} = undef;
	$self->{_time} = 0;
	
	
	return @re;
}


=head1 Accessor METHODS

=cut

=head2 fh

Get/set the file handle the progress bar is printed to.

=cut

sub fh{
	my ($self, $fh) = @_;
	$self->{fh} = $fh if defined($fh);
	return $self->{fh}
}

=head2 size

Get/set size.

=cut

sub size{
	my ($self, $size) = @_;
	$self->{size} = $size if defined($size);
	return $self->{size}
}

=head1 private METHODS

=cut

=head2 eta

=cut




=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
