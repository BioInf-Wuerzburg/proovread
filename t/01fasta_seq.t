#!/usr/bin/env perl

# $Id: Verbose.pm 644 2013-01-18 16:04:44Z dumps $

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";


#--------------------------------------------------------------------------#
=head2 sample data

=cut

my $head = <<HEAD;
>gi|129295|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)
HEAD

my $seq = <<SEQ;
QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAE
KMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTS
VLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHP
FLFLIKHNPTNTIVYFGRYWSP
SEQ

my $seq_raw = $head.$seq;

my $seq_obj_empty = {
	'seq_head' => '',
	'desc' => '',
	'byte_offset' => undef,
	'id' => '',
	'seq' => ''
};

my $seq_obj = {
	'seq_head' => '>gi|129295|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)',
	'desc' => 'GENE X PROTEIN (OVALBUMIN-RELATED)',
	'byte_offset' => undef,
	'id' => 'gi|129295|sp|P01013|OVAX_CHICK',
	'seq' => 'QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAEKMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTSVLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHPFLFLIKHNPTNTIVYFGRYWSP'
};

my $seq_string = <<SEQSTRING;
>gi|129295|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)
QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAEKMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTSVLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHPFLFLIKHNPTNTIVYFGRYWSP
SEQSTRING

#--------------------------------------------------------------------------#
=head2 load module

=cut

BEGIN { use_ok('Fasta::Seq'); }

my $class = 'Fasta::Seq';

#--------------------------------------------------------------------------#
=head2 new

=cut

# NOTE: cannot use is_deeply as it crashes on overloaded objects...
# empty
my $fa_empty = new_ok($class);

subtest 'new empty object' => sub{
	foreach my $attr(keys %$seq_obj_empty){
		is($fa_empty->{$attr}, $seq_obj_empty->{$attr}, "attribute $attr")
	}
};

# string
my $fa = Fasta::Seq->new($seq_raw);

subtest 'new from string' => sub{
	foreach my $attr(keys %$seq_obj){
		is($fa->{$attr}, $seq_obj->{$attr}, "attribute $attr")
	}
};

# hash
my $fa_hash = Fasta::Seq->new(
	seq => $seq,
	seq_head => $head,
);

subtest 'new from hash' => sub{
	foreach my $attr(keys %$seq_obj){
		is($fa_hash->{$attr}, $seq_obj->{$attr}, "attribute $attr")
	}
};

# clone
my $fa_clone = $fa->new();

subtest 'new clone' => sub{
	foreach my $attr(keys %$seq_obj){
		is($fa_clone->{$attr}, $seq_obj->{$attr}, "attribute $attr")
	}
};


#--------------------------------------------------------------------------#
=head2 Accessors

=cut

subtest 'generics: can and get' => sub{
	can_ok($class, 'seq_head');
	is($fa->seq_head(), $seq_obj->{seq_head}, 'seq_head() get');
	can_ok($class, 'id');
	is($fa->id(), $seq_obj->{id}, 'id() get');
	can_ok($class, 'desc');
	is($fa->desc(), $seq_obj->{desc}, 'desc() get');
	can_ok($class, 'seq');
	is($fa->seq(), $seq_obj->{seq}, 'seq() get');
	can_ok($class, 'byte_offset');
	is($fa->byte_offset(), $seq_obj->{byte_offset}, 'byte_offset() get');
};

subtest 'generics: set and dependencies' => sub{	
	# seq_head set and dependencies
	$fa->seq_head('>some new header');
	is($fa->{seq_head}, '>some new header', 'seq_head() set');
	$fa->seq_head('some new header');
	is($fa->{seq_head}, '>some new header', 'seq_head() set missing >');
	
	is($fa->id(), 'some', 'seq_head() set dependency id');
	is($fa->desc(), 'new header', 'seq_head() set dependency desc');
	
	# id set and dependencies
	$fa->id('foo');
	is($fa->id(), 'foo', 'id() set');
	is($fa->seq_head(), '>foo new header', 'id() set dependency seq_head');
	is($fa->desc(), 'new header', 'id() set dependency desc');
	
	# desc set and dependencies
	$fa->desc('bar');
	is($fa->desc(), 'bar', 'desc() set');
	is($fa->seq_head(), '>foo bar', 'desc() set dependency seq_head');
	is($fa->id(), 'foo', 'desc() set dependency id');

	$fa->desc('');
	is($fa->desc(), '', 'desc() reset');
	is($fa->seq_head(), '>foo', 'desc() reset dependency seq_head');
	is($fa->id(), 'foo', 'desc() reset dependency id');

	$fa->seq("AACT\nAGC");
	is($fa->seq(), 'AACTAGC', 'seq() set with \n trimming');

	$fa->byte_offset(5);
	is($fa->byte_offset(), 5, 'byte_offset() set');

};

subtest '$obj->desc_append' => sub{
	can_ok($class, 'desc_append');
	is($fa->desc_append('bla'), 'bla', 'desc_append() return value on empty desc');
	is($fa->desc(), 'bla', 'desc_append() on empty desc');
	is($fa->desc_append('blub'), 'bla blub', 'desc_append() return value on non-empty desc');
	is($fa->desc(), 'bla blub', 'desc_append() on non-empty desc');
};

subtest '$obj->string' => sub{
	can_ok($class, 'string');
	is($fa_clone->string, $seq_string, 'string() with one line sequence');
	is($fa_clone->string(70), $seq_raw, 'string() with linebreaks');
	is("$fa_clone", $seq_string, 'string() overload ""');
};

subtest '$obj->complement' => sub{
	can_ok($class, 'complement');
	$fa->seq('ATGCatgcAATTT');
	is($fa->complement(), 'TACGtacgTTAAA', 'complement()');
};

subtest '$obj->cat' => sub{
	can_ok($class, 'cat');
	my $seq1 = "ATGGGCA";
	my $seq2 = "ataNNaatgg";
	my $seq_app = $seq1.$seq2;
	my $seq_pre = $seq2.$seq1;
	my $fa1 = Fasta::Seq->new(">id1\n".$seq1);
	my $fa2 = Fasta::Seq->new(">id2 descr\n".$seq2);
	
	# Class Method
	is(Fasta::Seq->cat($fa1,$fa2)->seq, $seq_app, 'cat() class method, append obj');
	is(Fasta::Seq->cat($fa1,$fa2,1)->seq, $seq_pre, 'cat() class method, prepend obj');
	is(Fasta::Seq->cat($fa1,$seq2)->seq, $seq_app, 'cat() class method, append seq');
	is(Fasta::Seq->cat($seq2, $fa1)->seq, $seq_pre, 'cat() class method, prepend seq');
	
	# Object Method
	is($fa1->cat($fa2)->seq, $seq_app, 'cat() object method, append obj');
	is($fa1->cat($fa2, 1)->seq, $seq_pre , 'cat() object method, prepend obj');
	is($fa1->cat($seq2)->seq, $seq_app, 'cat() object method, append seq');
	is($fa1->cat($seq2, 1)->seq, $seq_pre , 'cat() object method, prepend seq');
	
	# Overload Method
	ok(overload::Method($fa1, '.'),'cat() overloads "."');
	my $fa_app_obj = $fa1.$fa2;
	is($fa_app_obj->seq, $seq_app, 'cat() overload ".", append obj');
	my $fa_app_seq = $fa1.$fa2;
	is($fa_app_seq->seq, $seq_app, 'cat() overload ".", append seq');
	my $fa3=$fa1->new;
	$fa3.=$fa2;
	is($fa3->seq, $seq_app, 'cat() overload ".=", append obj');
	$fa1.=$seq2;
	is($fa1->seq, $seq_app, 'cat() overload ".=", append seq');
	
};

subtest '$obj->reverse_complement' => sub{
	can_ok($class, 'base_content');
	$fa->seq('ATGCatgcAATTT');
	is($fa->reverse_complement(), 'AAATTgcatGCAT', 'reverse_complement()');
};


subtest '$obj->base_content' => sub{
	can_ok($class, 'base_content');
	$fa->seq('AATGGGCCGGNAAGGaatggc');
	is($fa->base_content('A'), 4, 'base_content(), single char');
	is($fa->base_content('GG'), 7, 'base_content(), mulitple chars');
};

subtest '$obj->substr_seq' => sub{
	can_ok($class, 'substr_seq');
	$fa->seq('ATGCatgcAATTT');
	is($fa->substr_seq(1,5)->seq, 'TGCat', 'substr_seq() offset, length');
	is($fa->substr_seq(1)->seq, 'TGCatgcAATTT', 'substr_seq() offset');
	is($fa->substr_seq(-4)->seq, 'ATTT', 'substr_seq() -offset');
	is($fa->substr_seq(1,-3)->seq, 'TGCatgcAA', 'substr_seq() offset, -length');
	is($fa->substr_seq(1,-1, 'XXX')->seq, 'AXXXT', 'substr_seq() offset, -length, replacement');
	
	my @subseqs = $fa->substr_seq([1,5], [1,-1, 'XXX']);
	is($subseqs[0]->seq, 'TGCat', 'substr_seq() batch offset, length');
	is($subseqs[1]->seq, 'AXXXT', 'substr_seq() batch offset, -length, replacement');
};

done_testing();

__END__


