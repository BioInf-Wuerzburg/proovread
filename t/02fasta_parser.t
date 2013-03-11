#!/usr/bin/env perl

# $Id: Verbose.pm 644 2013-01-18 16:04:44Z dumps $

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

BEGIN { use_ok('Fasta::Parser'); }

my $class = 'Fasta::Parser';
my ($fp, $fpr, $fpr_string);
my $fpr_dat_file = 'fasta.fa';

my $seq1_string = <<'FASTA';
>gi|129295|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)
QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAEKMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTSVLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHPFLFLIKHNPTNTIVYFGRYWSP
FASTA
my $seq2_string = <<SEQSTRING;
>gi|444439576|ref|NR_074891.1| Escherichia coli O157:H7 str. Sakai strain Sakai 16S ribosomal RNA, complete sequence
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGAGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTA
SEQSTRING
my $seq3_string = <<SEQSTRING;
>gi|436904082|gb|JX952207.1| Homo sapiens precursor microRNA HepG2-8, complete sequence
AGGGTCCAGCAGGGAGGGGTGGCAGGGATTCCCGCCTGTTACAGAGACTCCCCCACACCTTTCTCCCAACAGGTGTTCCCG
SEQSTRING


# new
subtest 'new and DESTROY' => sub{
	$fp = new_ok($class);
	
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
	$fpr = new_ok($class, [
		mode => '<', 
		file => $fpr_dat_file,
	]);
	#  default attributes
	is($fpr->{mode}, '<', 'Custom attribute "mode"');
	is($fpr->{file}, $fpr_dat_file, 'Custom attribute "file"');
	
		# large file
	# create large file
	my $fasta_large = ($seq1_string.$seq2_string.$seq3_string)x10000;
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

	is($fpr->next_seq->string, $seq1_string, 'Get $obj->next_seq()');	
	is($fpr->next_seq->string, $seq2_string, 'Get $obj->next_seq()');
	is($fpr->next_seq->string, $seq3_string, 'Get $obj->next_seq()');
	is($fpr->next_seq, undef, 'Get $obj->next_seq() eof');
	is($fpr->next_seq->string, $seq1_string, 'Get $obj->next_seq() restart after eof');
};

# check_format

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
	my $fpw_tmp_file = 'fasta.tmp';
	my $fpw = new_ok($class, [
		mode => '+>', 
		file => $fpw_tmp_file,
	]);
	
	#  default attributes
	is($fpw->{mode}, '+>', 'Writing mode "mode"');
	is($fpw->{file}, $fpw_tmp_file, 'Writer mode "file"');
	is(-e $fpw_tmp_file, 1, 'Create custom file'); # copy of STDIN needs to be greater than 3
	
	unlink $fpw_tmp_file;
};

# append_seq

# append_tell





done_testing();
