#!/usr/bin/env perl

# $Id: 03fastq_seq.t 59 2013-05-16 12:57:10Z s187512 $

use strict;
use warnings;

use Test::More;
use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";


#--------------------------------------------------------------------------#
=head2 load module

=cut

BEGIN { use_ok('Fastq::Seq'); }

my $Class = 'Fastq::Seq';

#--------------------------------------------------------------------------#
=head2 sample data

=cut


# create data file names from name of this <file>.t
(my $Dat_file = $FindBin::RealScript) =~ s/t$/dat/; # data
(my $Dmp_file = $FindBin::RealScript) =~ s/t$/dmp/; # data structure dumped
(my $Tmp_file = $FindBin::RealScript) =~ s/t$/tmp/; # data structure dumped

# slurp <file>.dat
my $Dat = do { local $/; local @ARGV = $Dat_file; <> }; # slurp data to string
my %Dat;
@Dat{qw(pb1 il1 crp1 crp2)} = split(/(?<=\n)(?=@)/, $Dat);

# eval <file>.dump
my %Dmp;
@Dmp{qw(empty pb1 il1 crp1 crp2)} = do "$Dmp_file"; # read and eval the dumped structure


#--------------------------------------------------------------------------#
=head1 ClassMethods

=cut

=head2 CheckFormat

=cut

subtest 'CheckFormat' => sub {
	can_ok($Class, 'CheckFormat');
	is($Fastq::Seq::CheckFormat, 0, 'Default $Fastq::Seq::CheckFormat');
	is(Fastq::Seq->CheckFormat(), 0, 'Get Fastq::Seq->CheckFormat()');
	is(Fastq::Seq->CheckFormat(1), 1, 'Set Fastq::Seq->CheckFormat()');
	is($Fastq::Seq::CheckFormat, 1, 'Check $Fastq::Seq::CheckFormat');
	is(Fastq::Seq->CheckFormat(0), 0, 'Reset Fastq::Seq->CheckFormat()');
	is($Fastq::Seq::CheckFormat, 0, 'Check $Fastq::Seq::CheckFormat');
};

#--------------------------------------------------------------------------#
=head2 new

=cut

# NOTE: cannot use is_deeply as it crashes on overloaded objects...
# empty
my ($fq_empty, $fq_il1, $fq_pb1, $fq_crp1);

$fq_empty = new_ok($Class);
subtest 'new empty object' => sub{
	foreach my $attr(keys %{$Dmp{empty}}){
		is($fq_empty->{$attr}, $Dmp{empty}{$attr}, "attribute $attr")
	}
};



# string
$fq_il1 = $Class->new($Dat{il1});
$fq_pb1 = $Class->new($Dat{pb1});

subtest 'new from string' => sub{
	foreach my $attr(keys %{$Dmp{il1}}){
		is($fq_il1->{$attr}, $Dmp{il1}{$attr}, "attribute $attr (il)");
		is($fq_pb1->{$attr}, $Dmp{pb1}{$attr}, "attribute $attr (pb)");
	}
	eval{ $Class->new(undef) };
	like($@, qr/^Fastq::Seq::new/, '$obj->new fatal on undef STRING');
	eval{ $Class->new("some\nstring\n") };
	like($@, qr/^Fastq::Seq::new/, '$obj->new fatal on corrupt STRING');
	eval{ $Class->new("\@id", "SEQU") };
	like($@, qr/^Fastq::Seq::new/, '$obj->new fatal on corrupt ARRAY');
	
};


# array
my $fq_il1_hash = $Class->new(	@{$Dmp{il1}}{qw(seq_head seq qual_head qual)} );
my $fq_pb1_hash = $Class->new(	@{$Dmp{pb1}}{qw(seq_head seq qual_head qual)} );

subtest 'new from array' => sub{
	foreach my $attr(keys %{$Dmp{il1}}){
		is($fq_il1_hash->{$attr}, $Dmp{il1}{$attr}, "attribute $attr (il)");
		is($fq_pb1_hash->{$attr}, $Dmp{pb1}{$attr}, "attribute $attr (pb)");
	}
};

# clone
my $fq_il1_clone = $fq_il1->new();
my $fq_pb1_clone = $fq_pb1->new();

subtest 'new clone' => sub{
	foreach my $attr(keys %{$Dmp{il1}}){
		is($fq_il1_clone->{$attr}, $Dmp{il1}{$attr}, "attribute $attr (il)");
		is($fq_pb1_clone->{$attr}, $Dmp{pb1}{$attr}, "attribute $attr (pb)");
	}
};


subtest 'check_format / new with CheckFormat(1)' => sub{
	
	eval{ $Class->new($Dat{il1})->check_format };
	unlike($@, qr/^Fastq::Seq::check_format/, '$obj->check_format non-fatal on valid string');
	eval{ $Class->new($Dat{crp1})->check_format };
	like($@, qr/^Fastq::Seq::check_format/, '$obj->check_format fatal on missing qual_head');
	eval{ $Class->new($Dat{crp2})->check_format };
	like($@, qr/^Fastq::Seq::check_format/, '$obj->check_format fatal on differing length in qual in seq');
	
	eval{ $Class->new( @{$Dmp{pb1}}{qw(seq_head seq qual_head qual)} )->check_format };
	unlike($@, qr/^Fastq::Seq::check_format/, '$obj->check_format non-fatal on valid array');
	eval{ $Class->new( @{$Dmp{pb1}}{qw(seq_head seq qual_head qual)}, phred_offset => 64 )->check_format };
	like($@, qr/^Fastq::Seq::check_format/, '$obj->check_format fatal on invalid phred_offset');
	
	$Class->CheckFormat(1);
	eval{ $Class->new($Dat{crp1}) }; # capture fatal
	like($@, qr/^Fastq::Seq::check_format/, 'new after $class->CheckFormat(1) fatal on corrupted data');
};


#--------------------------------------------------------------------------#
=head2 Accessors

=cut

subtest 'generics' => sub{
	can_ok($Class, 'seq_head');
	is($fq_il1->seq_head(), $Dmp{il1}->{seq_head}, 'seq_head() get');
	$fq_il1->seq_head('Bar/1 FOO');
	is($fq_il1->seq_head(), '@Bar/1 FOO', 'seq_head() set');
	
	can_ok($Class, 'seq');
	is($fq_il1->seq(), $Dmp{il1}->{seq}, 'seq() get');
	$fq_il1->seq("ATGCGACC\nNNATGA");
	is($fq_il1->seq(), 'ATGCGACCNNATGA', 'seq() set');
	
	can_ok($Class, 'qual_head');
	is($fq_il1->qual_head(), $Dmp{il1}->{qual_head}, 'qual_head() get');
	$fq_il1->qual_head("hmlgr");
	is($fq_il1->qual_head(), "+hmlgr", 'qual_head() set');
	
	can_ok($Class, 'qual');
	is($fq_il1->qual(), $Dmp{il1}->{qual}, 'qual() get');
	$fq_il1->qual("IIIII\nIGGGHHH");
	is($fq_il1->qual(), "IIIIIIGGGHHH", 'qual() set');
	
	can_ok($Class, 'phred_offset');
	is($fq_il1->phred_offset(), $Dmp{il1}->{phred_offset}, 'phred_offset() get');
	
	can_ok($Class, 'id');
	is($fq_il1->id(), 'Bar/1', 'id() get');
	$fq_il1->id('@Foo.2:0/1');
	is($fq_il1->seq_head(), '@Foo.2:0/1 FOO', 'id() set dependency seq_head');
	is($fq_il1->id(), 'Foo.2:0/1', 'id() set');
	
	can_ok($Class, 'id_no_picard');
	$fq_il1->seq_head('Bar/1 FOO');
	is($fq_il1->id_no_picard(), 'Bar', 'id_no_picard() get');
	$fq_il1->id_no_picard('Foo');
	is($fq_il1->id(), 'Foo/1', 'id_no_picard() set dependency seq_head');
	is($fq_il1->id_no_picard(), 'Foo', 'id_no_picard() set');
	$fq_il1->seq_head('Foo/Bar');
	is($fq_il1->id_no_picard(), 'Foo/Bar', 'id_no_picard() get');
	$fq_il1->id('Foo/Bar/1');
	is($fq_il1->id(), 'Foo/Bar/1', 'id_no_picard() set dependency seq_head');
	is($fq_il1->id_no_picard(), 'Foo/Bar', 'id_no_picard() set');
	$fq_il1->id_no_picard('Foo/Bar/1234');
	is($fq_il1->id(), 'Foo/Bar/1234/1', 'id_no_picard() set dependency seq_head');
	is($fq_il1->id_no_picard(), 'Foo/Bar/1234', 'id_no_picard() set');
	
	can_ok($Class, 'picard');
	$fq_il1->seq_head('Bar/1 FOO');
	is($fq_il1->picard(), '1', 'picard() get');

	can_ok($Class, 'desc');
	is($fq_il1->desc(), 'FOO', 'desc() get');
	$fq_il1->desc('Hi There');
	is($fq_il1->seq_head(), '@Bar/1 Hi There', 'desc() set dependency seq_head');
	is($fq_il1->desc(), 'Hi There', 'desc() set');
	$fq_il1->desc('');
	is($fq_il1->seq_head(), '@Bar/1', 'desc() reset dependency seq_head');
	is($fq_il1->desc(), '', 'desc() reset');

};

subtest '$obj->desc_append' => sub{
	$fq_il1->seq_head('@Foo/1');
	can_ok($Class, 'desc_append');
	is($fq_il1->desc_append('bla'), 'bla', 'desc_append() return value on empty desc');
	is($fq_il1->desc(), 'bla', 'desc_append() on empty desc');
	is($fq_il1->desc_append('blub'), 'bla blub', 'desc_append() return value on non-empty desc');
	is($fq_il1->desc(), 'bla blub', 'desc_append() on non-empty desc');
};


subtest '$obj->string' => sub{
	can_ok($Class, 'string');
	is($fq_pb1->string, $Dat{pb1}, 'string()');
	is("$fq_pb1", $Dat{pb1}, 'string() overload ""');
};


subtest '$obj->cat' => sub{
	can_ok($Class, 'cat');
	my $seq1 = "ATGGGCA";
	my $qual1 ="IIIIHHH";
	my $seq2 = "ataNNaatgg";
	my $qual2 ='ii$$II***_';
	my $seq_app = $seq1.$seq2;
	my $qual_app = $qual1.$qual2;
	my $seq_pre = $seq2.$seq1;
	my $qual_pre = $qual2.$qual1;
	
	my $fq1 = $Class->new("\@id1\n".$seq1."\n+\n".$qual1);
	my $fq2 = $Class->new("\@id2 descr\n".$seq2."\n+\n".$qual2);
	
	# Class Method
	my $fq_cat = $Class->cat($fq1,$fq2);
	is($fq_cat->seq, $seq_app, 'cat() class method, append obj seq');
	is($fq_cat->qual, $qual_app, 'cat() class method, append obj qual');
	
	$fq_cat = $Class->cat($fq1,$fq2,1);
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
	can_ok($Class, 'base_content');
	$fq_il1->seq('AATGGGCCGGNAAGGaatggc');
	is($fq_il1->base_content('A'), 4, 'base_content(), single char');
	is($fq_il1->base_content('GG'), 7, 'base_content(), mulitple chars');
};

subtest '$obj->substr_seq' => sub{
	can_ok($Class, 'substr_seq');
	$fq_il1->seq('ATGCatgcAATTT');
	$fq_il1->qual('IIGGHAGHIIIII');
	is($fq_il1->substr_seq(1,5)->seq, 'TGCat', 'substr_seq() offset, length');
	is($fq_il1->substr_seq(1,5)->qual, 'IGGHA', 'substr_seq() offset, length');

	is($fq_il1->substr_seq(1)->seq, 'TGCatgcAATTT', 'substr_seq() offset');
	is($fq_il1->substr_seq(1)->qual, 'IGGHAGHIIIII', 'substr_seq() offset');

	is($fq_il1->substr_seq(-4)->seq, 'ATTT', 'substr_seq() -offset');
	is($fq_il1->substr_seq(-4)->qual, 'IIII', 'substr_seq() -offset');

	is($fq_il1->substr_seq(1,-3)->seq, 'TGCatgcAA', 'substr_seq() offset, -length');
	is($fq_il1->substr_seq(1,-3)->qual, 'IGGHAGHII', 'substr_seq() offset, -length');

	is($fq_il1->substr_seq(1,-1, 'XXX', 'BBB')->seq, 'AXXXT', 'substr_seq() offset, -length, replacement');
	is($fq_il1->substr_seq(1,-1, 'XXX', 'BBB')->qual, 'IBBBI', 'substr_seq() offset, -length, replacement');
	
	my @subseqs = $fq_il1->substr_seq([1,5], [1,-1, 'XXX', 'BBB']);
	is($subseqs[0]->seq, 'TGCat', 'substr_seq() batch offset, length');
	is($subseqs[0]->qual, 'IGGHA', 'substr_seq() batch offset, length');

	is($subseqs[1]->seq, 'AXXXT', 'substr_seq() batch offset, -length, replacement');
	is($subseqs[1]->qual, 'IBBBI', 'substr_seq() batch offset, -length, replacement');

};

done_testing();

__END__


