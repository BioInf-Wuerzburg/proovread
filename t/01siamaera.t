#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

##-----------------------------------------------------------------------------##

=head2 sample data

=cut

# create data file names from name of this <file>.t
my $bin = $RealBin.'/../bin/siamaera'; # siaemaera.pl
(my $Fa_file = $FindBin::RealScript) =~ s/t$/fa/; # data
(my $Dat_file = $FindBin::RealScript) =~ s/t$/dat/; # data structure dumped
(my $Dmp_file = $FindBin::RealScript) =~ s/t$/dmp/; # data structure dumped
(my $Tmp_file = $FindBin::RealScript) =~ s/t$/tmp/; # data structure dumped
(my $pre = $FindBin::RealScript) =~ s/.t$//;

my ($Dat, %Dat, %Dmp);

if(-e $Dat_file){
	# slurp <file>.dat
	$Dat = do { local $/; local @ARGV = $Dat_file; <> }; # slurp data to string
	# %Dat = split("??", $Dat);
}

if(-e $Dmp_file){
    # eval <file>.dump
    %Dmp = do "$Dmp_file"; # read and eval the dumped structure
}


# Test cases
##----------------------------------------------------------------------------##
# head                      m140218_000044_42149_c100624412550000001823118308061424_s1_p0/3248/0_4880
# tail (eq to 1.)           m140218_000044_42149_c100624412550000001823118308061424_s1_p0/3248/4928_9784
# unsplit J                 m140218_000044_42149_c100624412550000001823118308061424_s1_p0/4501/10696_16013.1
# J > 300                   m140218_000044_42149_c100624412550000001823118308061424_s1_p0/5772/0_6624
# low qcov                  m140218_000044_42149_c100624412550000001823118308061424_s1_p0/6308/12181_16058
# ??                        m140218_000044_42149_c100624412550000001823118308061424_s1_p0/9724/0_4368
# low idy with xgap > 50    m140218_000044_42149_c100624412550000001823118308061424_s1_p0/6308/12181_16058
# unsplit J, low idy        m140218_000044_42149_c100624412550000001823118308061424_s1_p0/13299/13695_15712
# 4+ HSPs
# 3+ HSPs
##----------------------------------------------------------------------------##

like(qx($bin --version), qr/\d+\.\d+/, "exec");

print qx($bin --debug <$Fa_file 1>/dev/null);

done_testing()
