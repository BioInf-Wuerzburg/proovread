#!/usr/bin/env perl

# $Id: 02fasta_parser.t 55 2013-05-15 11:41:39Z s187512 $

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

BEGIN { use_ok('Fasta::Parser'); }

my $Class = 'Fasta::Parser';
my ($fp, $fpr, $fpr_string);

#--------------------------------------------------------------------------#

=head2 sample data

=cut

(my $Dat_file = $FindBin::RealScript) =~ s/t$/dat/; # data
(my $Dmp_file = $FindBin::RealScript) =~ s/t$/dmp/; # data structure dumped
(my $Tmp_file = $FindBin::RealScript) =~ s/t$/tmp/; # data structure dumped

my $Dat = do { local $/; local @ARGV = $Dat_file; <> }; # slurp data to string
my @Dat = split(/(?<=\n)(?=>)/, $Dat);

my @Dmp = do "$Dmp_file"; # read and eval the dumped structure


# new
subtest 'new and DESTROY' => sub{
	$fp = new_ok($Class);
	
	#  default attributes
	is($fp->{mode}, '<', 'Default attribute "mode"');
	is($fp->{file}, undef, 'Default attribute "file"');
	is($fp->{_buffer}, undef, 'Default attribute "_buffer"');
	ok(fileno $fp->{fh} > 2, 'Default attribute "fh"'); # copy of STDIN needs to be greater than 3
	
	# DESTROY closes fh
	my $fh = $fp->{fh};
	$fp = undef;
	is(fileno($fh), undef, 'Autoclose filehandle');
	# read from file
	$fpr = new_ok($Class, [
		mode => '<', 
		file => $Dat_file,
	]);
	#  default attributes
	is($fpr->{mode}, '<', 'Custom attribute "mode"');
	is($fpr->{file}, $Dat_file, 'Custom attribute "file"');
	
		# large file
	# create large file
	my $fasta_large = ($Dat)x10000;
	$fpr_string = Fasta::Parser->new(file => \$fasta_large);
	
};


subtest 'filehandle' => sub{
	# fh
	can_ok($fpr, 'fh');
	can_ok($fpr, 'is_fh');
	my $fprh = $fpr->{fh};
	is($fpr->fh, $fprh, 'Get $obj->fh()');
	is($fpr->is_fh, 0, '$obj->is_fh() FILE');
	ok($fpr->is_fh('FILE'), '$obj->is_fh("FILE")');
	
	is($fpr->fh(\*STDIN), \*STDIN, 'Set $obj->fh()');
	is($fpr->is_fh, 1, '$obj->is_fh() PIPE');
	ok($fpr->is_fh('PIPE'), '$obj->is_fh("PIPE")');
	
	is($fpr->fh($fprh), $fprh , 'Reset $obj->fh() FILE');
	is($fpr->is_fh, 0, '$obj->is_fh()');

	is($fpr_string->is_fh(), 2, '$obj->is_fh() SCALAR');	
	ok($fpr_string->is_fh("SCALAR"),'$obj->is_fh("SCALAR")');	
	
};

# next_seq
subtest 'next_seq' => sub{
	can_ok($fpr, 'next_seq');

	for(my $i=0; $i<@Dat; $i++){
		is($fpr->next_seq->string, $Dat[$i], 'Get $obj->next_seq()');	
	}
	is($fpr->next_seq, undef, 'Get $obj->next_seq() eof');
	is($fpr->next_seq->string, $Dat[0], 'Get $obj->next_seq() restart after eof');
};

# check_format
subtest 'check_format' => sub{
	can_ok($fpr, 'check_format');
	is(Fasta::Parser->new(file => \("Not a fasta record"))->check_format, 
		undef, '$obj->check_format on not FASTA is undef');
	my $fp = Fasta::Parser->new(file => $Dat_file);
	isa_ok($fp->check_format, $Class, '$obj->check_format on FASTA');
	isa_ok($fp->check_format, $Class, '$obj->check_format on FASTA');
};

# seek

# sample_seqs
subtest 'sample_seqs' => sub{
	can_ok($fpr, 'sample_seqs');
	my @sample_seqs = $fpr->sample_seqs(2);
	is(scalar @sample_seqs, 2, '$obj->sample_seqs()');
	isa_ok($_, 'Fasta::Seq') for @sample_seqs;
	@sample_seqs = $fpr->sample_seqs(3);
	is(scalar @sample_seqs, 3, '$obj->sample_seqs()');
	isa_ok($_, 'Fasta::Seq') for @sample_seqs;
	@sample_seqs = $fpr->sample_seqs(10);
	is(scalar @sample_seqs, 3, '$obj->sample_seqs()');
	isa_ok($_, 'Fasta::Seq') for @sample_seqs;
	
	@sample_seqs = $fpr_string->sample_seqs(10);
	is(scalar @sample_seqs, 10, '$obj->sample_seqs() large file');
	isa_ok($_, 'Fasta::Seq') for @sample_seqs;
};

# guess_seq_length

# guess_seq_count

# check_fh_is_pipe


subtest 'new write parser' => sub{
	# writer
	my $fpw = new_ok($Class, [
		mode => '+>', 
		file => $Tmp_file,
	]);
	
	#  default attributes
	is($fpw->{mode}, '+>', 'Writing mode "mode"');
	is($fpw->{file}, $Tmp_file, 'Writer mode "file"');
	is(-e $Tmp_file, 1, 'Create custom file'); # copy of STDIN needs to be greater than 3
	
	unlink $Tmp_file;
};

# append_seq

# append_tell





done_testing();
