#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;

my $cfg_core_file = "../proovread.cfg";
my %cfg = do $cfg_core_file;

my @TASKS = @{$cfg{'mode-tasks'}{'pacbio-pre'}};

print Dumper(cfg('long-reads'));
print Dumper(cfg('coverage'));
print Dumper(cfg('sr-pre-fraction'));
print Dumper(cfg('sr-pre-fraction', 'shrimp-pre-1') );
print Dumper(cfg('sr-pre-fraction', 1));
print Dumper(cfg('sr-pre-fraction', -2));
print Dumper(scalar cfg('mode-tasks', 'pacbio-pre') );


# HASH
# Problem: there always needs to be the default (0) or a specific task to 
# retrieve stuff

sub cfg{
	my ($k, $t) = @_;
	return undef unless exists $cfg{$k};
	return $cfg{$k} unless ref $cfg{$k};
	
	if(ref $cfg{$k} eq "ARRAY"){
		return wantarray ?	@{$cfg{$k}} : $cfg{$k};
	}
	
	if(! exists $cfg{$k}{DEF}){
		return wantarray ?	%{$cfg{$k}} : $cfg{$k};
	}
	
	my %p = %{$cfg{$k}};
	return undef unless %p;
	
	my $v = exists $p{DEF} ? $p{DEF} : undef;
	
	if(defined $t){
		if($t =~ /[^-0-9]/){
			$v = $p{$t} if exists $p{$t};
		}else{
			$v = $p{$TASKS[$t]} if exists $p{$TASKS[$t]};
		}
	}
	return $v;
}

# ARRAY
# Pro: Default can be first entry of %2 array, rest goes to task hash
# Con: No normal ARRAYs can be used anymore

sub cfg_array{
	my ($k, $t) = @_;
	return undef unless exists $cfg{$k};
	return $cfg{$k} if (! ref $cfg{$k} || ref $cfg{$k} eq 'HASH');
	
	# assume ARRAY
	# first entry, if %2 is default
	my @p = @{$cfg{$k}};
	my $def = shift @p if @p%2; # any global
	
	return $def unless @p; # only global
	
	my %p = @p;
	return undef unless %p;
	
	my $v = defined $def ? $def : undef; # in case the task has no settings and there is no global
	
	if(defined $t){
		if($t =~ /[^-0-9]/){
			$v = $p{$t} if exists $p{$t};
		}else{
			$v = $p{$TASKS[$t]} if exists $p{$TASKS[$t]};
		}
	}
	return $v;
}


