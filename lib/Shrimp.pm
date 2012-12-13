package Shrimp;

use warnings;
use strict;

# $Id$

use File::Temp;

# preference libs in same folder over @INC
use lib './';

use Verbose;

our $VERSION = '0.05';
our ($REVISION) = '$Revision$' =~ /(\d+)/;
our ($MODIFIED) = '$Date$' =~ /Date: (\S+\s\S+)/;


$|++;



##------------------------------------------------------------------------##

=head1 NAME 

Shrimp.pm

=head1 DESCRIPTION

SHRiMP (2.2.0) gmapper interface.

=head1 SYNOPSIS

  use Shrimp;

=cut

=head1 CHANGELOG

=head2 0.05

=over

=item [Change] Preference libs in same folder over @INC

=item [Change] Added svn:keywords

=back

=over 12

=item 0.04 [Thomas Hackl 2012-10-29]

Takes now an 'out' file argument to write data to. If specified, 
 the output can not be read on the fly, but a filehandle to read
 the complete result is provided, after the run has finished.

=item 0.03 [Thomas Hackl 2012-10-26]

Added C<< $shrimp->status >> method to retrieve the current status message.

=item 0.02

Removed IPC::Open3. Does not work with gmapper, do know why. STDOUT
 is now captured via pipe, STDERR via temporary file

=item 0.01

Initial module. Provides Constructor and generic accessor and run methods.

=back

=cut

=head1 TODO

=over

=item Tests

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
	format => "[{TIME_SHORT}] [Shrimp] {MESSAGE}\n"
);

our $VB = Verbose->new(
	line_delim => "\\\n",
	line_width => 80
);

##------------------------------------------------------------------------##

=head1 Constructor METHOD


=head2 new

  out => 'FILE'    # file argument to write data to. If specified, 
                   #  the output can not be read on the fly and the run cannot be 
                   #  canceled manually or by timeout. Still a filehandle to read
                   #  the complete result is provided, after the run has finished.
  log => 'FILE'    # A file to write the log to, defaults to 'shrimpXXXXX.log', 
                   #  with X being random numbers
  verbose => BOOL  # Default TRUE
  timeout => INT   # Number of seconds before canceling run by timeout
  bin => 'gmapper-ls'

=cut


sub new {
	my ($class) = shift;
	my $self = {
		# defaults
		command => undef,
		bin => 'gmapper-ls',
		verbose => 1,
		timeout => 0,
		out => undef,
		log => undef,
		# overwrites
		@_,
		
		# protected privates
	};

	$V->{report_level} = $self->{verbose};
    # use defaults of called subClass 
    bless $self, $class;	
    
    $self->_build_command();
    
    return $self;
}




sub DESTROY{
	my $self = shift;
	close $self->{_error_reader} if $self->{_error_reader};
}


##------------------------------------------------------------------------##

=head1 Public METHODS

=head2 run

Public Method. Runs the blast search calling L<_run_ipc_open> or L<_run_pipe_open> if C<no_ipc = 1> is set. 
Returns the blast object.

  my $shrimp->run;

  --------------------------------------------------------------------------------
  gmapper: LETTER SPACE (454,Illumina/Solexa,etc.).
  SHRiMP 2.2.3
  [GCC 4.6.3; CXXFLAGS="-O3 -DNDEBUG -mmmx -msse -msse2 -fopenmp -Wall -Wno-deprecated -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -DGIT_VERSION=missing_git_version"]
  --------------------------------------------------------------------------------
  usage: gmapper-ls [options/parameters] { <r> | -1 <r1> -2 <r2> } <g1> <g2>...
     <r>                  Reads filename, paired or unpaired
     <r1>                 Upstream reads filename
     <r2>                 Downstream reads filename
     <g1> <g2>...         Space seperated list of genome filenames
  Parameters:
     -s/--seeds           Spaced Seed(s)                (default: 11110111101111,
                                                         1111011100100001111,
                                                         1111000011001101111)
     -o/--report          Maximum Hits per Read         (default: 10)
        --max-alignments  Max. align. per read  (0=all) (default: 0)
     -w/--match-window    Match Window Length           (default: 140.00%)
     -n/--cmw-mode        Match Mode                    (default: unpaired:2 paired:4)
     -l/--cmw-overlap     Match Window Overlap Length   (default: 90.00%)
     -a/--anchor-width    Anchor Width Limiting Full SW (default: 8; disable: -1)

     -S/--save            Save Genome Proj. in File     (default: no)
     -L/--load            Load Genome Proj. from File   (default: no)
     -z/--cutoff          Projection List Cut-off Len.  (default: 4294967295)

     -m/--match           SW Match Score                (default: 10)
     -i/--mismatch        SW Mismatch Score             (default: -15)
     -g/--open-r          SW Gap Open Score (Reference) (default: -33)
     -q/--open-q          SW Gap Open Score (Query)     (default: -33)
     -e/--ext-r           SW Gap Extend Score(Reference)(default: -7)
     -f/--ext-q           SW Gap Extend Score (Query)   (default: -3)
     -r/--cmw-threshold   Window Generation Threshold   (default: 55.00%)
     -h/--full-threshold  SW Full Hit Threshold         (default: 50.00%)

     -N/--threads         Number of Threads             (default: 1)
     -K/--thread-chunk    Thread Chunk Size             (default: 1000)

     -p/--pair-mode       Paired Mode                   (default: none)
     -I/--isize           Min and Max Insert Size       (default: 0,1000)
        --longest-read    Maximum read length           (default: 1000)
     -1/--upstream        Upstream read pair file
     -2/--downstream      Downstream read pair file
        --un              Dump unaligned reads to file
        --al              Dump aligned reads to file
        --read-group      Attach SAM Read Group name
        --sam-header      Use file as SAM header
        --single-best-mapping Report only the best mapping(s), this is not strata (see README)
        --all-contigs     Report a maximum of 1 mapping for each read.
        --no-mapping-qualities Do not compute mapping qualities
        --insert-size-dist Specifies the mean and stddev of the insert sizes
        --no-improper-mappings (see README)
        --trim-front      Trim front of reads by this amount
        --trim-end        Trim end of reads by this amount
        --trim-first      Trim only first read in pair
        --trim-second     Trim only second read in pair
        --min-avg-qv      The minimum average quality value of a read
        --progress        Display a progress line each <value> reads. (default 100000)




=cut


sub run {
	my ($self) = @_;
	
	$self->_run_pipe_open();
	
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
 started shrimp run. Checks for errors, removes tempfiles.

  my $shrimp->finish;

=cut

sub finish{
	my ($self) = @_; 
	unless (ref $self || ref $self !~ /^Shrimp/ ){
		die "Shrimp not initialized!\n";
	}
	
	my $err;
	if($self->{_tmp_error_log} && -e $self->{_tmp_error_log}){
		open(my $err, '<', $self->{_tmp_error_log});
		
		$self->{_error} = join('', <$err>);
		close $err;
		
		# redirect error reader
		open(my $err2, '<', \($self->{_error}));
		$self->{_error_reader} = $err2;
		
		# delete stderr file 
		unlink $self->{_tmp_error_log} unless $self->{'log'};
	}
	
	if($?){
		$self->{_status} = 'canceled';
	}else{
		# set status to 
		$self->{_status} = 'done';
	}

	$V->verbose($self->{_status});
	
	return $self;
}



=head2 cancel

Public Method. Cancels a running shrimp run based on a passed query process id 
or the internally stored process id of the shrimp object;

  my $shrimp->cancel(<message>);
  
  my Shrimp::Shrimp->cancel(pid, <message>);

=cut

sub cancel {
	my ($pid, $msg);
	# object method
	if(ref (my $me = shift)){
		$pid = $me ->{_pid};
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

	$V->verbose($msg || "Shrimp run canceled");
}


##------------------------------------------------------------------------##

=head1 Accessor METHODS

=cut

=head2 command

Get/Set the command.

=cut

sub command{
	my ($self, $command) = @_;
	$self->{command} = $command if $command;
	return $self->{command};
}

=head2 log

Get/Set the log file name.

=cut

sub log{
	my ($self, $log) = @_;
	$self->{'log'} = $log if $log;
	return $self->{'log'};
}


=head2 log_string

Get the log of the run as STRING.

=cut

sub log_string{
	my ($self) = @_;
	return $self->{_error};
}

=head2 oh

Get the output handle to read the output while generated by shrimp

=cut

sub oh{
	my ($self) = @_;
	return $self->{_result_reader};
}

=head2 logh

Get the log handle to read the log generated by shrimp

=cut

sub logh{
	my ($self) = @_;
	return $self->{_error_reader};
}

=head2 status

Get the status message ("running", "done", "canceled") of the shrimp run.

=cut

sub status{
	my ($self) = @_;
	return $self->{_status};
}


##------------------------------------------------------------------------##

=head1 Private Methods

=cut

=head1 _build_command

=cut

sub _build_command{
	my $self = shift;
	my @command;
	
	# binary
	if(-e $self->{bin}){
		push @command, $self->{bin}
	}else{
		my $which = "which ".$self->{bin}." 2>/dev/null";
		my $bin = qx($which);
		if($? >> 8){ # 0: success, 1: no bin, 2: option error
			$V->exit('Binary not found or not executable ('.$self->{bin}.')')
		}
		$bin =~ s/\n//g;
		if(-e $bin){
			$self->{bin} = $bin;
			push @command, $bin;
		}else{
			$V->exit('Binary not found or not executable ('.$self->{bin}.')')
		}
	}
	
	# params
	while(my($k, $v) = each %$self){
		next unless $k =~ /^-/;
		next if grep{$k eq $_}qw(-ref -1 -2); 
		# flag only is undef or '', NOT '0' !!!
		push @command, (defined($v) && $v ne '') ? ($k, $v) : $k;
	}
	
	# reads
	$V->exit("reads missing (-1, -2)") unless $self->{'-1'};
	if($self->{'-2'}){
		push @command, '-1', $self->{'-1'}, '-2', $self->{'-2'};
	}else{
		push @command, $self->{'-1'}
	}
	
	# genome
	$V->exit("Reference missing (-ref)") unless $self->{'-ref'};
	push @command, ref $self->{'-ref'} ? join(" ", @{$self->{'-ref'}}) : $self->{'-ref'};
	
	$self->{command} = join(" ", @command);	
}




=head2 _run_pipe_open

Private Method. Dispatches the blast call with a C<$pid = open my $rdr, '-|', "cmd"> construct to the shell 
and stores the result reader as well as a pid to cancel the run (actually not the pid returned by the open cmd) in the blast object.
Set C<no_ipc = 1> to prevent default L</_run_ipc_open> method. 

=cut

sub _run_pipe_open {
	my ($self) = shift;

	# logging
	my $tmp_error_log = $self->{'log'};
	unless($tmp_error_log){
		$tmp_error_log = File::Temp->new(
			TEMPLATE => "shrimpXXXXXX", 
			SUFFIX => ".log", 
			UNLINK => 0
		)->filename;
	};
	
	
	# This 'son of a bitch' construct obviously starts two processes.
	# Looks like one for the shell command with parent_pid +1
	# and one thats listening to the STOUT pipe coming back, with parent_pid +2 - ?, but as far as I know closer than +10.
	# killing the first process doesnt stop the blast, might even make the listener wait forever
	# killing the listener seams to be sufficient to end both processes, sth like waitpid in first, is parent
	
	my $rdr;
	my $open_pid;
	if($self->{out}){
		system($self->{command}.' 1> '.$self->{out}.' 2> '.$tmp_error_log);
		open ($rdr, '<', $self->{out}) or $V->exit($!); 
	
	}else{
		$open_pid = open ($rdr, $self->{command}.' 2> '.$tmp_error_log.' |') 
			or $V->exit($!);

		$self->{_status} = 'running';
		$V->verbose($self->{_status}.", pid: $open_pid");
	
		# little hack to verify the listeners pid
		# get the entry for the all running gmapper processes
		# search for entry, where open pid is parent and then take the childs id as new id
		#sleep(1);
		my $ps = qx(ps -eF | grep -v grep | grep $self->{bin});
		($self->{_pid}) = $ps  =~ /^\S+\s+(\d+)\s+$open_pid/m;
	}

	$self->{_result_reader} = $rdr;

	open(my $err, '<', $tmp_error_log) or $V->exit($!);
	$self->{_error_reader} = $err;
	$self->{_tmp_error_log} = $tmp_error_log;
	
	return $self;
}

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


