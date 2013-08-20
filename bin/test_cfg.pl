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
print Dumper(cfg('sr-pre-fraction', 'shrimp-pre-1'));
print Dumper(cfg('sr-pre-fraction', 1));
print Dumper(cfg('sr-pre-fraction', -2));



sub cfg{
	my ($k, $t) = @_;
	return undef unless exists $cfg{$k};
	return $cfg{$k} unless ref $cfg{$k};
	
	if(ref $cfg{$k} eq "ARRAY"){
		return @{$cfg{$k}};
	}
	
	my %p = %{$cfg{$k}};
	return undef unless %p;
	
	my $v = exists $p{0} ? $p{0} : undef;
	
	if(defined $t){
		if($t =~ /[^-0-9]/){
			$v = $p{$t} if exists $p{$t};
		}else{
			$v = $p{$TASKS[$t]} if exists $p{$TASKS[$t]};
		}
	}
	return $v;
}
