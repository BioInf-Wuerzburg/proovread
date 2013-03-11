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

my $pb_raw = <<'PB';
@m120727_181124_42164_c100331452550000001523019509161202_s2_p0/47/117_1776
CACCGGAAACAAAGCGGACTTTGCTTCTTTTCCATCCATTTTCATTGCGTTTTC
+
(&,)..(#&.&&&'!'"&"++"/.+-$,'%&%-%.+/&/*(+/&,+.(./)+,+
PB

my($pb_sh, $pb_seq, $pb_qh, $pb_qual) = split("\n", $pb_raw);

my $il_raw = <<'IL';
@Foo.2:0/1 bar boo
CAGGCTAAAACAAACAGGGCCTAGGTTTTTGGTCATCGGGATAAAATGACGGTT
+
AAAAAVYYYYYYYYYYYYYbbbbbbbbb_________bbbbYYYYYYYYYYbbb
IL

my($il_sh, $il_seq, $il_qh, $il_qual) = split("\n", $il_raw);


my $seq_obj_empty = {
	'seq_head' => '@',
	'seq' => '',
	'qual_head' => '+',
	'qual' => '',
	'phred_offset' => undef,
};

# 1
my $seq_obj_il = {
	'seq_head' => $il_sh,
	'qual_head' => $il_qh,
	'phred_offset' => undef,
	'seq' => $il_seq,
	'qual' => $il_qual
};

# 2
my $seq_obj_pb = {
	'seq_head' => $pb_sh,
	'qual_head' => $pb_qh,
	'phred_offset' => undef,
	'seq' => $pb_seq,
	'qual' => $pb_qual
};

#--------------------------------------------------------------------------#
=head2 load module

=cut

BEGIN { use_ok('Fastq::Seq'); }

my $class = 'Fastq::Seq';

#--------------------------------------------------------------------------#
=head2 new

=cut

# NOTE: cannot use is_deeply as it crashes on overloaded objects...
# empty
my $fq_empty = new_ok($class);
subtest 'new empty object' => sub{
	foreach my $attr(keys %$seq_obj_empty){
		is($fq_empty->{$attr}, $seq_obj_empty->{$attr}, "attribute $attr")
	}
};




# string
my $il = $class->new($il_raw);
my $pb = $class->new($pb_raw);

subtest 'new from string' => sub{
	foreach my $attr(keys %$seq_obj_il){
		is($il->{$attr}, $seq_obj_il->{$attr}, "attribute $attr (il)");
		is($pb->{$attr}, $seq_obj_pb->{$attr}, "attribute $attr (pb)");
	}
};


# hash
my $il_hash = $class->new(	$il_sh, $il_seq, $il_qh, $il_qual );
my $pb_hash = $class->new(	$pb_sh, $pb_seq, $pb_qh, $pb_qual );


subtest 'new from array' => sub{
	foreach my $attr(keys %$seq_obj_il){
		is($il_hash->{$attr}, $seq_obj_il->{$attr}, "attribute $attr (il)");
		is($pb_hash->{$attr}, $seq_obj_pb->{$attr}, "attribute $attr (pb)");
	}
};


# clone
my $il_clone = $il->new();
my $pb_clone = $pb->new();

subtest 'new clone' => sub{
	foreach my $attr(keys %$seq_obj_il){
		is($il_clone->{$attr}, $seq_obj_il->{$attr}, "attribute $attr (il)");
		is($pb_clone->{$attr}, $seq_obj_pb->{$attr}, "attribute $attr (pb)");
	}
};

#--------------------------------------------------------------------------#
=head2 Accessors

=cut

subtest 'generics' => sub{
	can_ok($class, 'seq_head');
	is($il->seq_head(), $seq_obj_il->{seq_head}, 'seq_head() get');
	$il->seq_head('Bar/1 FOO');
	is($il->seq_head(), '@Bar/1 FOO', 'seq_head() set');
	
	can_ok($class, 'seq');
	is($il->seq(), $seq_obj_il->{seq}, 'seq() get');
	$il->seq("ATGCGACC\nNNATGA");
	is($il->seq(), 'ATGCGACCNNATGA', 'seq() set');
	
	can_ok($class, 'qual_head');
	is($il->qual_head(), $seq_obj_il->{qual_head}, 'qual_head() get');
	$il->qual_head("hmlgr");
	is($il->qual_head(), "+hmlgr", 'qual_head() set');
	
	can_ok($class, 'qual');
	is($il->qual(), $seq_obj_il->{qual}, 'qual() get');
	$il->qual("IIIII\nIGGGHHH");
	is($il->qual(), "IIIIIIGGGHHH", 'qual() set');
	
	can_ok($class, 'phred_offset');
	is($il->phred_offset(), $seq_obj_il->{phred_offset}, 'phred_offset() get');
	
	can_ok($class, 'id');
	is($il->id(), 'Bar/1', 'id() get');
	$il->id('@Foo.2:0/1');
	is($il->seq_head(), '@Foo.2:0/1 FOO', 'id() set dependency seq_head');
	is($il->id(), 'Foo.2:0/1', 'id() set');
	
	can_ok($class, 'id_no_picard');
	$il->seq_head('Bar/1 FOO');
	is($il->id_no_picard(), 'Bar', 'id_no_picard() get');
	$il->id_no_picard('Foo');
	is($il->id(), 'Foo/1', 'id_no_picard() set dependency seq_head');
	is($il->id_no_picard(), 'Foo', 'id_no_picard() set');
	$il->seq_head('Foo/Bar');
	is($il->id_no_picard(), 'Foo/Bar', 'id_no_picard() get');
	$il->id('Foo/Bar/1');
	is($il->id(), 'Foo/Bar/1', 'id_no_picard() set dependency seq_head');
	is($il->id_no_picard(), 'Foo/Bar', 'id_no_picard() set');
	$il->id_no_picard('Foo/Bar/1234');
	is($il->id(), 'Foo/Bar/1234/1', 'id_no_picard() set dependency seq_head');
	is($il->id_no_picard(), 'Foo/Bar/1234', 'id_no_picard() set');
	
	can_ok($class, 'picard');
	$il->seq_head('Bar/1 FOO');
	is($il->picard(), '1', 'picard() get');

	can_ok($class, 'desc');
	is($il->desc(), 'FOO', 'desc() get');
	$il->desc('Hi There');
	is($il->seq_head(), '@Bar/1 Hi There', 'desc() set dependency seq_head');
	is($il->desc(), 'Hi There', 'desc() set');
	$il->desc('');
	is($il->seq_head(), '@Bar/1', 'desc() reset dependency seq_head');
	is($il->desc(), '', 'desc() reset');

};


subtest '$obj->desc_append' => sub{
	$il->seq_head('@Foo/1');
	can_ok($class, 'desc_append');
	is($il->desc_append('bla'), 'bla', 'desc_append() return value on empty desc');
	is($il->desc(), 'bla', 'desc_append() on empty desc');
	is($il->desc_append('blub'), 'bla blub', 'desc_append() return value on non-empty desc');
	is($il->desc(), 'bla blub', 'desc_append() on non-empty desc');
};


subtest '$obj->string' => sub{
	can_ok($class, 'string');
	is($pb->string, $pb_raw, 'string()');
	is("$pb", $pb_raw, 'string() overload ""');
};


subtest '$obj->cat' => sub{
	can_ok($class, 'cat');
	my $seq1 = "ATGGGCA";
	my $qual1 ="IIIIHHH";
	my $seq2 = "ataNNaatgg";
	my $qual2 ='ii$$II***_';
	my $seq_app = $seq1.$seq2;
	my $qual_app = $qual1.$qual2;
	my $seq_pre = $seq2.$seq1;
	my $qual_pre = $qual2.$qual1;
	
	my $fq1 = $class->new("\@id1\n".$seq1."\n+\n".$qual1);
	my $fq2 = $class->new("\@id2 descr\n".$seq2."\n+\n".$qual2);
	
	# Class Method
	my $fq_cat = $class->cat($fq1,$fq2);
	is($fq_cat->seq, $seq_app, 'cat() class method, append obj seq');
	is($fq_cat->qual, $qual_app, 'cat() class method, append obj qual');
	
	$fq_cat = $class->cat($fq1,$fq2,1);
	is($fq_cat->seq, $seq_pre, 'cat() class method, prepend obj seq');
	is($fq_cat->qual, $qual_pre, 'cat() class method, append obj qual');
	
	# Object Method
	is($fq1->cat($fq2)->seq, $seq_app, 'cat() object method, append obj seq');
	is($fq1->cat($fq2)->qual, $qual_app, 'cat() object method, append obj qual');
	is($fq1->cat($fq2, 1)->seq, $seq_pre , 'cat() object method, prepend obj seq');
	is($fq1->cat($fq2, 1)->qual, $qual_pre , 'cat() object method, prepend obj qual');
	
	# Overload Method
	ok(overload::Method($fq1, '.'),'cat() overloads "."');
	$fq_cat = $fq1.$fq2;
	is($fq_cat->seq, $seq_app, 'cat() overload ".", append obj seq');
	is($fq_cat->qual, $qual_app, 'cat() overload ".", append obj qual');
	
	my $fq3=$fq1->new;
	$fq3.=$fq2;
	is($fq3->seq, $seq_app, 'cat() overload ".=", append obj');
	
};



subtest '$obj->base_content' => sub{
	can_ok($class, 'base_content');
	$il->seq('AATGGGCCGGNAAGGaatggc');
	is($il->base_content('A'), 4, 'base_content(), single char');
	is($il->base_content('GG'), 7, 'base_content(), mulitple chars');
};

subtest '$obj->substr_seq' => sub{
	can_ok($class, 'substr_seq');
	$il->seq('ATGCatgcAATTT');
	$il->qual('IIGGHAGHIIIII');
	is($il->substr_seq(1,5)->seq, 'TGCat', 'substr_seq() offset, length');
	is($il->substr_seq(1,5)->qual, 'IGGHA', 'substr_seq() offset, length');

	is($il->substr_seq(1)->seq, 'TGCatgcAATTT', 'substr_seq() offset');
	is($il->substr_seq(1)->qual, 'IGGHAGHIIIII', 'substr_seq() offset');

	is($il->substr_seq(-4)->seq, 'ATTT', 'substr_seq() -offset');
	is($il->substr_seq(-4)->qual, 'IIII', 'substr_seq() -offset');

	is($il->substr_seq(1,-3)->seq, 'TGCatgcAA', 'substr_seq() offset, -length');
	is($il->substr_seq(1,-3)->qual, 'IGGHAGHII', 'substr_seq() offset, -length');

	is($il->substr_seq(1,-1, 'XXX', 'BBB')->seq, 'AXXXT', 'substr_seq() offset, -length, replacement');
	is($il->substr_seq(1,-1, 'XXX', 'BBB')->qual, 'IBBBI', 'substr_seq() offset, -length, replacement');
	
	my @subseqs = $il->substr_seq([1,5], [1,-1, 'XXX', 'BBB']);
	is($subseqs[0]->seq, 'TGCat', 'substr_seq() batch offset, length');
	is($subseqs[0]->qual, 'IGGHA', 'substr_seq() batch offset, length');

	is($subseqs[1]->seq, 'AXXXT', 'substr_seq() batch offset, -length, replacement');
	is($subseqs[1]->qual, 'IBBBI', 'substr_seq() batch offset, -length, replacement');

};

done_testing();

__END__


