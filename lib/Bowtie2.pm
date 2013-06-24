package Bowtie2;

use warnings;
use strict;

# $Id: Bowtie2.pm 55 2013-05-15 11:41:39Z s187512 $

use File::Temp;
use Data::Dumper;

use IPC::Open3;

# preference libs in same folder over @INC
use lib './';

use Verbose;

our $VERSION = '0.01';


$|++;



##------------------------------------------------------------------------##

=head1 NAME 

Bowtie2.pm

=head1 DESCRIPTION

Bowtie2 interface.

=head1 SYNOPSIS

  use Bowtie2;
  
  my $bowtie2 = Bowtie2->new(
    ref => 'genome.fa'
    reads => 'reads.fa'
  );
  
  my $bowtie2->bowtie2;
  
  # read output on the fly
  use Sam::Parser;
  my $sp = Sam::Parser->new(
    fh => $bowtie2->oh
  );
  
  while(my $aln = $sp->next_aln()){
    # do something with the output
  }
  
  $bowtie2->finish;


=cut

=head1 CHANGELOG

head2 0.01

=item [INIT] Provides Constructor and generic accessor and run methods.

=back

=cut



##------------------------------------------------------------------------##

=head1 Class ATTRIBUTES

=cut

=head2 $V

Verbose messages are handled using the Verbose.pm module. To 
 customize verbose message behaviour, overwrite the attribute with
 another Verbose object created with the Verbose module.

=cut

our $V = Verbose->new(
	format => "[{TIME_SHORT}] [Bowtie2] {MESSAGE}\n"
);

our $VB = Verbose->new(
	line_delim => "\\\n",
	line_width => 80
);

=head1 Class METHODS

=cut

=head2 Param_join (HASHREF, JOIN=STRING)

Joins a HASHREF to a parameter string with JOIN [" "], ignoring keys with 
 undef values and creating flag only values for ''.

  Bowtie2->Param_join(HASHREF); # join with space
  Bowtie2->Param_join(HASHREF, join => "\n") # join with newline

=cut

sub Param_join{
	my $proto = shift;
	my $params = shift @_;
	my $p = {
		'join' => " ",
		'ignore' => [],
		@_
	};
	
	my %ignore;
	@ignore{@{$p->{ignore}}} = (1)x scalar @{$p->{ignore}}
		if @{$p->{ignore}};
	# params
	my @params;
	my $paramstring;
	foreach my $k (sort keys %$params){
		next if exists $ignore{$k};
		my $v = $params->{$k};
		# flag only is '', NOT '0' !!!
		next unless defined ($v);
		push @params, ($v ne '') ? ($k, $v) : $k;
	}
	$paramstring .= join($p->{'join'}, @params);
	return $paramstring;
}




##------------------------------------------------------------------------##

=head1 Constructor METHOD


=head2 new

  out => 'FILE'    # file argument to write data to. If specified, 
                   #  the output can not be read on the fly and the run cannot be 
                   #  canceled manually or by timeout. Still a filehandle to read
                   #  the complete result is provided, after the run has finished.
  log => 'FILE'    # A file to write the log to, defaults to 'bowtie2XXXXX.log', 
                   #  with X being random numbers
  verbose => BOOL  # Default TRUE
  timeout => INT   # Number of seconds before canceling run by timeout
  bin => 'bowtie2'
                   # use explicit path, in case bowtie2 cannot be found in
                   #  the path, use bowtie2 -C for color space
  ref => [FILE, FILE,...]
                   # file or list of files containing reference sequences
  reads => 'FILE'  # file containing short reads
  mates => 'FILE'  # file containing mate reads in case paired mapping is performed




=cut


sub new {
	my ($class) = shift;
    
	my $self = {
		# defaults
		bowtie2_bin => 'bowtie2',
		bowtie2_opt => {},
		bowtie2_build_bin => 'bowtie2-build',
		bowtie2_build_opt => {},
		path => '',
		verbose => 1,
		timeout => 0,
		'ref' => undef,
		'pre' => undef,
		'-x' => undef,
		'-1' => undef,
		'-2' => undef,
		'-U' => undef,
		'-S' => undef,
		# overwrites
		@_,
		# "protected"
		_status => 'initialized',
		_pid => undef,
		_timeout_pid => undef,
		_stdin => undef,
		_stdout => undef,
		_stderr => undef,
		# protected privates
	};

	$V->{report_level} = $self->{verbose};
    # use defaults of called subClass 
    bless $self, $class;	
    
    return $self;
}




sub DESTROY{
	my $self = shift;
}


##------------------------------------------------------------------------##

=head1 Public METHODS

=head2 run

Start bowtie2 run with parameters provided to new method. Returns bowtie2 object.

  my $bowtie2 = $bowtie2->run;

=cut

sub bowtie2_build {
	my $self = shift;
	
	my $p = {
		# "globals"
		'ref' => $self->{ref},
		'pre' => $self->{pre},
		# bowtie2_build specific
		%{$self->{bowtie2_build_opt}},
		# custom overwrite
		@_
	};
	
	# store the overwritten bowtie2_build opts - this looks somehow dangerous
	$self->{bowtie2_build_opt} = $p; 

	die("Current status '".$self->status()."'. 'finish' prior to any other new action")
		if $self->status =~ /^running/;

	my $opts = $self->opt2string("bowtie2_build", ignore=>[qw(ref pre)]);
	
	# let open3 do its magic :)	
	use Symbol 'gensym'; 
	$self->{_stderr} = gensym;
	$self->{_pid} = open3(
		$self->{_stdin},
		$self->{_stdout},
		$self->{_stderr},
		$self->path ? $self->path.'/'.$self->{bowtie2_build_bin} : $self->{bowtie2_build_bin},
		($opts ? $opts : ()),
		$p->{ref},
		$p->{pre}
	);

	$self->{_status} = "running bowtie2_build";

	# fork timeout process to monitor the run and cancel it, if necessary
	# child
	if(!$self->{out} && $self->{timeout} && !( $self->{_timeout_pid} = fork()) ){
		# fork error
		if ( not defined $self->{_timeout_pid} ){ die "couldn't fork: $!\n"; }
		$self->_timeout(); 
		# childs exits either after timeout canceled blast or blast run has finished
	}
	# parent - simply proceeds
	
	return $self
}

sub bowtie2 {
	my $self = shift;
		my $p = {
		# "globals"
		'-1' => $self->{'-1'},
		'-2' => $self->{'-2'},
		'-x' => $self->{'-x'},
		'-S' => $self->{'-S'},
		'-U' => $self->{'-U'},
		# bowtie2_build specific
		%{$self->{bowtie2_opt}},
		# custom overwrite
		@_
	};
	
	# store the overwritten bowtie2_build opts - this looks somehow dangerous
	$self->{bowtie2_opt} = $p; 

	die("Current status '".$self->status()."'. 'finish' prior to any other new action")
		if $self->status =~ /^running/;
		
	# let open3 do its magic :)	
	use Symbol 'gensym'; 
	$self->{_stderr} = gensym;
	$self->{_pid} = open3(
		$self->{_stdin},
		$self->{_stdout},
		$self->{_stderr},
		$self->path ? $self->path.'/'.$self->{bowtie2_bin} : $self->{bowtie2_bin},
		$self->opt2string("bowtie2")
	);
	
	$self->{_status} = "running bowtie2";
	
	# fork timeout process to monitor the run and cancel it, if necessary
	# child
	if(!$self->{out} && $self->{timeout} && !( $self->{_timeout_pid} = fork()) ){
		# fork error
		if ( not defined $self->{_timeout_pid} ){ die "couldn't fork: $!\n"; }
		$self->_timeout(); 
		# childs exits either after timeout canceled blast or blast run has finished
	}
	# parent - simply proceeds
	
	return $self
}


=head2 finish

Public Method. Waits for the finishing/canceling of a 
 started bowtie2 run. Checks for errors, removes tempfiles.

  my $bowtie2->finish;

=cut

sub finish{
	my ($self) = @_; 

	unless (ref $self || ref $self ne "Bowtie2" ){
		die "Bowtie2 not initialized!\n";
	}
	
	
	# to make sure, bowtie2 is finished, read its STDOUT until eof
	unless($self->{out}){
		my $tmp;
		1 while read($self->{_stdout},$tmp,10000000);
	}

	waitpid($self->{_pid}, 0);
	
	if($?){
		unless($self->status =~ /^canceled/){ # set by cancel method
			my $exitval = $? >> 8;
			$self->{_status} =~ s/running/exited($exitval)/;
		}
	}else{
		# set status to 
		$self->{_status} =~ s/running/finished/;
	}

	return $self;
}



=head2 cancel

Public Method. Cancels a running bowtie2 run based on a passed query process id 
or the internally stored process id of the bowtie2 object;

  my $bowtie2->cancel(<message>);
  
  Bowtie2->cancel(pid, <message>);

=cut

sub cancel {
	my ($pid, $msg);
	# object method
	if(ref (my $me = shift)){
		$pid = $me ->{_pid};
		$me->{_status} =~ s/^\w+/canceled/;

		$msg = shift;
	}
	# class method
	else{
		($pid, $msg) = @_;
	}
	
	unless( kill ('1', $pid) ){
		# TODO: cancel pid does not exist
#		$pid." doesnt exist -> probably already finished\n");
	};

}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 opt2string

Get a stringified version of the specified parameter for a command.

  $bowtie2_opt = $self->opt2string("bowtie2");
  $bowtie2_build_opt = $self->opt2string("bowtie2_build");
  

=cut

sub opt2string{
	my ($self, $cmd, @more) = @_;
	die("unknown command") unless exists $self->{$cmd.'_opt'};
	return Bowtie2->Param_join($self->{$cmd.'_opt'}, @more) || '';
}

=head2 path

Get/Set the path to the binaries.

=cut

sub path{
	my ($self, $path) = @_;
	$self->{path} = $path if defined($path);
	return $self->{path};
}

=head2 status

Get the current status of the bowtie object.

=cut

sub status{
	my ($self) = @_;
	return $self->{_status};
}

=head2 stderr

Get the filehandle to the stderr stream. Only available while status is 
'running'.

=cut

sub stderr{
	my ($self) = @_;
	return $self->{_stderr};
}

=head2 stdout

Get the filehandle to the stdout stream. Only available prior to C<< $bowtie2->finish >>.

=cut

sub stdout{
	my ($self) = @_;
	return $self->{_stdout};
}

=head2 log_string

Get the log of the run as STRING.

=cut

sub log_string{
	my ($self) = @_;
	return $self->{_error};
}

=head2 oh

Get the output handle to read the output while generated by bowtie2

=cut

sub oh{
	my ($self) = @_;
	return $self->{_result_reader};
}

=head2 logh

Get the log handle to read the log generated by bowtie2

=cut

sub logh{
	my ($self) = @_;
	return $self->{_error_reader};
}


##------------------------------------------------------------------------##

=head1 Private Methods

=cut

=head2 timeout

Private Method. Initialized by run if C<timeout> > than 0;

=cut

sub _timeout {
	my ($self) = @_; 
	my $time = 0;
	# set sleep to default 2 seconds if not specified
	if(not defined $self->{_sleep} ){ 
		# set sleep time to never be higher than timeout, but maximal 2 seconds
		$self->{_sleep} = $self->{timeout} < 2 ? $self->{timeout} : 2;
	}
	$V->verbose('timeout of '.$self->{timeout}.'s activated');
	while(my $pid = kill ('0', $self->{ _pid}) ){ 
		if($time > $self->{timeout}){
			$self->cancel( 'Canceled by timeout '.$self->{timeout}."s" );
			exit(0); 
		}else{
			$time += $self->{_sleep};	
			sleep($self->{_sleep});
		}
	}
	exit(0);
}




=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;


