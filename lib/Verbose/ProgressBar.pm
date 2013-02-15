package Verbose::ProgressBar;

use warnings;
use strict;

# $Id$

# preference libs in same folder over @INC
use lib '../';

use Verbose;

our $VERSION = '0.08';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;

$|++;



##------------------------------------------------------------------------##

=head1 NAME Progress

=head1 DESCRIPTION

Collection of generic function for progress calculation and visualisation

=cut

=head1 SYNOPSIS

=cut

=head1 TODO

=over

=item Docu

=item SYNOPSIS

=item unknown size

=item very small sizes (< bin)

=back

=cut

=head1 CHANGELOG

=head2 0.08

=over

=item [BugFix] TSS is now always printed in case finish is called.

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
	
	if (!defined ($value) && defined($self->{_fh})){
		$self->{value} = tell($self->{_fh});
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
			$eta = sprintf(" TTS %02d:%02d:%02d", (gmtime($time-$self->{_start_time}))[2,1,0])
		}elsif($eta = ($self->{_bins}-$bin) * $self->{_bin_time}){
			$eta = sprintf(" ETA %02d:%02d:%02d", (gmtime($eta))[2,1,0])
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
