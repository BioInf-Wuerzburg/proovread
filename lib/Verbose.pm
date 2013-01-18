package Verbose;

use warnings;
use strict;

# $Id$

# preference libs in same folder over @INC
use lib './';

our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;
our $VERSION = '0.07';



##------------------------------------------------------------------------##

=head1 NAME 

Verbose

=head1 DESCRIPTION

Collection of generic function to create verbose messages.

=cut

=head1 SYNOPSIS

  my $v = Verbose->new(
    format => "[{CALLER} line:{LINE} {TIME_ELAPSED}] {MESSAGE}\n",
  );

  $v->verbose('Hello World');
  $v->nline();
  $v->verbose('Hello World again');
    # [main line 13 00:00:00] Hello World
    # 
    # [main line 15 00:00:00] Hello World again

  my $v_bash = Verbose->new(
    report_level => 0,
    fh => \*STDOUT,
    line_delim => "\\\n",
    line_width => 50
  );

  $v_bash->verbose('foo bar');
    # prints nothing because level (default 1) is higher then report_level.

  $v_bash->hline(level => 0);
  $v_bash->verbose(
    message =>'perl /my/long/long/path/to/my/perl/script.pl --foo="this parameter" --bar="that parameter" --noverbose --prefix="some/thing"',
    level => 0
  );
    # #------------------------------------------------#
    # perl /my/long/long/path/to/my/perl/script.pl --foo\
    # ="this parameter" --bar="that parameter" --noverbo\
    # se --prefix="some/thing"
  
  package Some::Package 
  
  some_sub();
   
  # [Some::Package line:13 00:00:00] FATAL ERROR
  #   Package: Some::Package
  #   File:    C:/projects/biosrc/sandbox/bin/Some/Package.pm
  #   Line:    13
  #   Sub:     Some::Package::some_sub
  # and exits the running program with "exit 1".
  
  sub some_sub{
    $v->exit('FATAL ERROR')
  }

=cut

=head1 TODO

=cut

=head1 CHANGELOG

=head2 0.08

=over

=item [Feature] Added attribute C<exit_val> with default 255 to Verbose objects.

=back

=head2 0.07

=over

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=over

=item 0.06 [Thomas Hackl 2012-10-26]

Removed STDERR aliasing, causes problems, if different Verbose instances
 are used simultaniously.

=item 0.05 [Thomas Hackl 2012-10-24]

Minor POD corrections.

=item 0.04 [Thomas Hackl 2012-10-22]

Default fh changed from STDERR to ">&STDERR" (copy of STDERR) which
 can and will be safely closed if object leaves scope.

=item 0.03 [Thomas Hackl]

Added more detailed trace info. Adapted C<exit()> method to display this
 information by default.

=item 0.02

Added C<Humanize()> class method to create human readable numbers.

=item 0.01

Initial Verbose module. Provides C<new(), verbose(), exit(), 
 nline() and hline()> method, generic accessor methods and 
 formatting templates.

=back

=cut

=head1 Constructor METHOD

=head2 Attributes

=over 10

=item report_level

Maximum level for a message to be printed, default 1.

=item level

Default level of a message, if not set explicitly, default 1.

=item format

Default template, default "{MESSAGE}\n".

=item fh

Filehandle to print to, default \*STDERR

=item line_width

Number of characters after which line break is introduced.
Default 0 => unlimited.

=item line_delim

Line delimiting character(s). Default "\n".

=item start_time

Start time for the ELAPSED_TIME computation. Default is system time.

=back

=cut

sub new{
	my $class = shift;
	
	# defaults=
	my $self = {
		fh => \*STDERR,
		level => 1,
		report_level => 1,
		format => "[{TIME_SHORT}] {MESSAGE}\n",	# default template is class template
		format_exit => undef,
		start_time => time(),	# safe time() of creation for elapsed templates
		line_delim => "\n",
		line_width => undef,
		exit_val => 255,
		@_	# overwrite
	};
	
	if ($self->{file}){
		my $fh;
		open ( $fh , $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
		$self->{fh} = $fh;
	}
	
	$self->{format_exit} = "{HLINE}".$self->{format}."{TRACE}\n" unless $self->{format_exit};
	
	bless $self, $class;
}


=head1 Class ATTRIBUTES

=cut

=head2 Verbose::Templates

Collection of standart verbose templates. Use the names in '{}' within the 
format specification to be replaced with the resulting strint. 
Change/Create with Syntax 
  $Verbose::Templates::<TEMPLATE> = sub{ 
  	my $self = shift;
  	# do stuff; 
  	return STRING
  }.
  
$_[0] gives access to the object, the method is called on and thus can be 
its attributes can be used within the sub.

Predefined templates are:

=over 12

=item MESSAGE

The formatted message given to verbose.

=item CALLER

The name of the caller module, "main" if main program.

=item LINE

The line number from which the verbose message was sent.

=item TIME_FULL

Timestamp including date and day.

=item TIME_SHORT

Timestamp HH:MM:SS.

=item TIME_ELAPSED

Elapsed time since the verbose object has been created ~ run time of the
 program, HH:MM:SS.

=item HLINE

A horizontal line of format #----//--# of line_width or default 80 chars.

=back

=cut

our $Templates = {
	MESSAGE => sub{
		$_[0]->{message}
	},
	CALLER => sub{
		@{$_[0]->{trace_1}} ?  $_[0]->{trace_1}[0] : 'main'
	},
	LINE => sub{
		$_[0]->{trace_0}[2]
	},
	TIME_SHORT => sub{
		sprintf("%02d:%02d:%02d", (localtime)[2,1,0])
	},
	TIME_FULL => sub{
		sprintf("%s", scalar(localtime))
	},
	TIME_ELAPSED => sub{
		sprintf("%02d:%02d:%02d", (gmtime(time() - $_[0]->{start_time} ))[2,1,0])
	}, 
	HLINE => sub{
		return '#'.('-'x($_[0]->{line_width} ? $_[0]->{line_width} -2 : 78))."#\n" ;
	},
	TRACE => sub{
		my $self = shift;
		# package, filename, line, subroutine, hasargs
		sprintf("\tExited at '%s'\n\t'%s', line %s\n", @{$self->{trace_0}}[0,1,2]).
		(@{$self->{trace_1}} ? sprintf("\tLast call '%s'\n\t'%s', line %s", @{$self->{trace_1}}[3,1,2]) : '');
	}
};

=head2 $Verbose::SuffixTable

Suffix strings for 1000x substitution used by C<Verbose::Humanize()>

  12 => 'T',
  9 => 'G',
  6 => 'M',
  3 => 'k',
  0 => '',
  -3 => 'm',
  -6 => 'µ',
  -9 => 'n',
  -12 => 'p'


=cut

our $SuffixTable = {
	12 => 'T',
	9 => 'G',
	6 => 'M',
	3 => 'k',
	0 => ' ',
	-3 => 'm',
	-6 => 'µ',
	-9 => 'n',
	-12 => 'p'
};


=head1 Class METHODS

=cut

=head2 Humanize

Turn a number to a rounded, human readable string using 1000x suffixes.
 Takes a number, and an optional a precision parameter (default 2), returns 
 a rounded, formatted STRING. The maximum length of the returned STRING 
 computes as precision+4 (sign, precision+1 digits, floating ".", suffix
 string).

  print Verbose->Humanize(24534, precision=>2);
    # 2.45k

=cut

sub Humanize{
	my $proto = shift;
	my $class = ref $proto || $proto;
	
	my $p = {
		precision => 2,
		@_%2 ? (value => shift, @_) : @_
	};
	
	defined $p->{value} || return;
	
	
	my ($sig, $num, $exp) = 	sprintf("%.$p->{precision}e", $p->{value}) =~ m/
		^([+-]?)		# sign
		([\d.]+)		# num
		e([+-]\d+)		# exp
	/x;
	my $mod = $exp%3;

	if($mod){
		$num =~ tr/.//d;
		unless((my $pre = $mod+1) > $p->{precision}){
			$num =~s/(\d{$pre})/$1\./;
		};
	}

	join("", $sig, $num, $SuffixTable->{$exp-$mod});
}


=head1 Object METHODS

=cut

=head2 verbose

Print verbose message according to specified level.

=cut

sub verbose{
	my $self = shift;
	my $p = {
		%$self,
		trace_0 => [caller(0)],
		trace_1 => [caller(1)],
		@_%2 ? (message => shift, @_) : @_	# overwrite + odd number of params -> first is message
	}; 
	
	$p->{message} = '' unless defined ($p->{message});
	
	if ( $p->{level} <= $p->{report_level} ){
		
		# call format function with self and params, overwrite self
		(my $string = $p->{format}) =~ s/{(.+?)}/
			(exists ($Verbose::Templates->{$1}) && defined ($Verbose::Templates->{$1}))
				? &{$Verbose::Templates->{$1}}({%$self, %$p})
				: "{$1 - undefined replacement}";
		/gxe;
		
		# line width
		if(my $lw = $p->{line_width}){
			$lw -= length($p->{line_delim}) - 1;
			my $ndl = length($p->{line_delim}) || 0;
			$string = join ("\n", map{ 
				my $o=-$ndl; 
				substr($_, $o, 0, $p->{line_delim}) while (($o+=$lw+$ndl) < length($_)-1 );
				$_;
			}(split (/\n/, $string, -1)));
		}
		
		print {$self->{fh}} $string;
	}
}


=head2 hline

Print horizontal line.

=cut

sub hline{
	my $self = shift;
	
	my $p = {
		%$self,
		@_%2 ? (message => shift, @_) : @_	# overwrite + odd number of params -> first is message
	}; 
	
	if ( $p->{level} <= $p->{report_level} ){
		print {$self->{fh}} &{$Verbose::Templates->{HLINE}}($p); 
	}
}

=head2 nline

Print new line.

=cut

sub nline{
	my $self = shift;
	
	my $p = {
		%$self,
		@_%2 ? (message => shift, @_) : @_	# overwrite + odd number of params -> first is message
	}; 
	
	if ( $p->{level} <= $p->{report_level} ){
		print {$self->{fh}} "\n"; 
	}
}


=head2 exit

Print critical verbose message and exit with error (exit 1). Level is
 0, which means it is executed even if verbose messages are turned of.

=cut

sub exit{
	my $self = shift;
	$self->verbose(
		@_, 
		level => 0,
		trace_0 => [caller(0)],
		trace_1 => [caller(1)],
		format => $self->{format_exit}, # overwrite format to verbose temp.
	);
	exit $self->{exit_val};
}



=head1 Accessor METHODS

=head2 fh

Return/set the file handle the messages are directed to.

=cut

sub fh{
	my $self = shift;
	$self->{fh} = shift if @_ ;
	return $self->{fh};	
}

=head2 file

Return/set the file the messages are printed to.

=cut

sub file{
	my $self = shift;
	if (@_){
		$self->{file} = shift;
		my $fh;
		open ( $fh , $self->{file}) or die __PACKAGE__."file(): '".$self->{file}."', $!";
		$self->{fh} = $fh;
	}
	return $self->{file};	
}


=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut


1;