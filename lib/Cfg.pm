package Cfg;

=head1 NAME

Cfg.pm

=head1 DESCRIPTION

Read a simple config in perl hash style from file. Inspired from
 L<http://www.perlmonks.org/?node_id=464358>.

=head1 SYNOPSIS

  use Cfg;
  
  Cfg->Read_Cfg("path/to/my.cfg") || die "Couldn't read *.cfg";
  
  my $log_file = $Cfg::Cfg{log_file};

=head1 CHANGELOG

see git log.

=head1 TODO

=cut

#-----------------------------------------------------------------------------#

=head1 Include

=cut


use warnings;
no warnings 'qw';
use strict;

use Carp;
use Log::Log4perl qw(:easy);

use File::Spec;
use File::Copy;
use File::Basename;
use File::Which;

#use Data::Dumper;
#$Data::Dumper::Sortkeys = 1;



#-----------------------------------------------------------------------------#

=head1 Globals

=cut


my $L = Log::Log4perl::get_logger();


#-----------------------------------------------------------------------------#

=head2 Class Attributes

=cut

our $VERSION = 0.01;

our %Cfg;


#-----------------------------------------------------------------------------#

=head2 Aliases (for backward compatibility)

=cut

*Read_Cfg = \&Read;
*Copy_Cfg = \&Copy;


#-----------------------------------------------------------------------------#

=head1 Class Methods

=cut

=head2 Read

Read config from simple perl syntax config file (LIST context). Returns config
 , dies on error.

Simple *.cfg format:

  #--------------------------------#
  log => "path/to/somewhere",
  
  servers => [
    "foo", "bar"
  ],
  
  more => {
    complex => [qw( s t u f f )]
  },
  #--------------------------------#

=cut

sub Read{
	my ($class,$cfg_file, $subtree) = @_;

	unless(-f $cfg_file){
		$L->logdie("Cannot find config file '$cfg_file'");
	}

	%Cfg = do($cfg_file);
	
        if ($@) {
            $L->logdie("Failed to read config '$cfg_file' - $@");
        }

	if($subtree){
	    unless(exists $Cfg{$subtree}){
		$L->info('$Cfg{'.$subtree.'} does not exists, reading nothing');
		return;
	    }
	    $L->logdie('$Cfg{'.$subtree.'} is not a HASHREF - currently only hash subtrees are supported') unless ref $Cfg{$subtree} eq 'HASH';
	    
	    return %{$Cfg{$subtree}};
	}

	return %Cfg;
}


=head2 Copy

Create a copy of a config file at target. Default target is
<CWD>/basename(<SOURCE>).

  my $cfg_file_copy = Cfg->Copy_Cfg($cfg_file);
  my $cfg_file_copy = Cfg->Copy_Cfg($cfg_file, $new_cfg_file);


=cut

sub Copy{
	my ($class,$cfg_file, $new_cfg_file) = @_;

	unless(-f $cfg_file){
		$L->logdie("Cannot find config file '$cfg_file'");
	}

	$new_cfg_file = basename($cfg_file) unless $new_cfg_file;
	copy($cfg_file, $new_cfg_file) or $L->logdie("Creatring config failed: $!");
	
	return File::Spec->rel2abs($new_cfg_file);
}


=head2 Check_binaries

Test fatally whether a set of binaries are exported and/or existent
and executable.

  $Cfg->Check_binares(qw(/path/to/some/bin some-exported-bin ..)) or die "Missing binaries";

=cut

sub Check_binaries{
    my $class = shift;
    foreach my $bin (@_){
        unless(-e $bin && -x $bin){
            if(my $fbin = which($bin)){
                $L->logdie("Binary '$fbin' not executable") unless -e $fbin && -x $fbin;
            }else{
                $L->logdie("Binary '$bin' neither in PATH nor executable");
            }
        }
        
        $L->debug("Using binaries: ", which($bin) || $bin);
    }
    return 1;
}



=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut

1;


















