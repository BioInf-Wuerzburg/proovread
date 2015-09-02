package Sam::Phase;

use strict; 
use warnings;

use Sam::Seq;

use Data::Dumper;

our $VERSION = '0.1.0';


=head1 NAME

Sam::Phase.pm

=head1 DESCRIPTION

Perform variant phasing for Sam::Seq data.

=cut

=head1 SYNOPSIS

  use Sam::Seq;
  use Sam::Phase;

  my $sp = Sam::Parser->new(file => "index.bam");

  while (my $ss = $sp->next_seq){
    $ss->call_variants;
    $ss->phase_variants;
    ..
  }

=cut


=head1 Class ATTRIBUTES / METHODS

Get/Set ...

=cut

# Class attributes
our %ATTR_CLASS = (
    FragMinOvl => 2,
    FragMinLen => 3,
);

=head2 FragMinLen [3]

Minimum length of a haplotype fragments (number of variant positions) linked on
a read to be considered for phasing

=head2 FragMinOvl [2]

Minimum overlap of haplotype fragments from different reads to be considered for
joining.


##------------------------------------------------------------------------##

=head2 Sam::Seq::phase_variants

Greedy, seed based, heuristic approach for read-backed haplotype phasing. Based
on Levy et al, 2007. Well suited for phasing of long read data.

=cut

my $support_matrix;

sub Sam::Seq::phase_variants{
    my $self = shift;

    die (((caller 0)[3]).": ->call_variants() required\n") unless @{$self->{vars}};

    $self->{haplotype} = '';
        
    my $vars = $self->{vars};
    my @vpos = $self->variant_positions;

    return unless @vpos;
    my @stable_ref = $self->stable_ref_states;
    my @stable_ref_vars = @stable_ref[@vpos];
    
    my $vi_min = 0;

    my $c=0;
    my $frags = [];
    my $max_end;
    my $max_len;
    my $max_len_i;
    my $i;
    my @haplotypes;
    $support_matrix = {0 => [], 1 => []};

    foreach my $aln ( $self->alns(1) ) {

        my $o = $aln->pos -1;
        my @s = $aln->seq_states(gap => '');
        
        # first var pos index
        while ($vi_min < @vpos && $vpos[$vi_min] < $o) { $vi_min++ };
        last if $vi_min > $#vpos;        

        # last var pos indexq
        my $vi_max = $vi_min;

        while ($vi_max < @vpos && $vpos[$vi_max] < $o+@s) { $vi_max++ };
        $vi_max--;
        next if $vi_max - $vi_min <= Sam::Phase->FragMinLen;

        my @v_aln;
        my $hpl;
        
        foreach ( $vi_min .. $vi_max ) {
            my $v_aln_pos = $vpos[$_]-$o;
            my $v;
            my $l;

            if ( ($l = length($stable_ref_vars[$_])) == 1 ) {
                $v = $s[$v_aln_pos];
            }else {
                my $v_aln_to = $v_aln_pos+$l-1;
                $v_aln_to = $#s if $v_aln_to > $#s;
                $v = join("", @s[$v_aln_pos..$v_aln_to]);
            }

            push @v_aln, $v;
            $hpl.= $v eq $stable_ref_vars[$_] ? 1 : 0
        }

        my $frag = {
            seq => join ("", @v_aln),
            pos => $vi_min,
            end => $vi_max,
            hpl => $hpl,
        };
        
        #print $aln;
        #print " "x $vi_min, $frag->{hpl},"\n";
        #last if ++$c >= 100;

        if ( @$frags && $max_end - $frag->{pos} > Sam::Phase->FragMinOvl-2 ) { # >3bp overlap
            $i++;
            $max_end = $frag->{end} if $frag->{end} > $max_end;
            if ( length($frag->{seq}) > $max_len ){
                $max_len = length($frag->{seq});
                $max_len_i = $i;
            }
            push @$frags, $frag;
            next;
        }
        
        if (@$frags) {
            push @haplotypes, $self->compute_local_haplotype(frags => $frags, seed_idx => $max_len_i);
    }
        
        $i = 0;
        $frags = [$frag];
        $max_end = $frag->{end};
        $max_len = length($frag->{seq});
        $max_len_i = $i;
        
    }

    
    # in case we exited for loop on short frag
    if (@$frags) {
        push @haplotypes, $self->compute_local_haplotype(frags => $frags, seed_idx => $max_len_i);
    }

    # merge into global haplotype
    my $global_haplotype = '-' x @vpos;
    foreach ( @haplotypes ) {
        substr($global_haplotype, $_->{pos}, length($_->{hpl}), $_->{hpl});
    }

    #print $global_haplotype,"\n";

    # refine by majority vote of supporters
    my $refine = '';
    #my $_0 = '';
    #my $_1 = '';
    for ( my $i=0; $i<@vpos; $i++ ) {
        my $s0 = $support_matrix->{0}[$i] // 0;
        my $s1 = $support_matrix->{1}[$i] // 0;
        my $h = substr($global_haplotype, $i, 1);

        if ($h) {
            if ( $s1 < $s0 ) {
                $refine.= '*';
                substr($global_haplotype, $i, 1, "0");
            } else {
                $refine.= ' ';
            }
        } else {
            if ( $s0 < $s1 ) {
                $refine .= '*';
                substr($global_haplotype, $i, 1, "1");
            } else {
                $refine .= ' ';
            }
        }
    }

    $self->{haplotype} = $global_haplotype; # TODO: needs to have 1>>0.

    # flip vars according to ref
    foreach (@vpos) {
        my $ref = $stable_ref[$_];
        my @v = @{$vars->[$_]};

        my $i;
        for ($i=0;$i<@v;$i++) {
            last if $ref eq $v[$i];
        }

        # ref doesn't match primary but existing var
        if ($i && $i<@v) { # flip var order
            @{$self->{vars}[$_]}[0,$i] = @{$self->{vars}[$_]}[$i,0];
            @{$self->{freqs}[$_]}[0,$i] = @{$self->{vars}[$_]}[$i,0];
            @{$self->{probs}[$_]}[0,$i] = @{$self->{probs}[$_]}[$i,0];
        }
    }
}


sub Sam::Seq::haplotype_variants{
    my $self = shift;

    die (((caller 0)[3]).": ->phase_variants() required\n") unless defined ($self->{haplotype});

    my @flip_var;
    my @hpl = split(//, $self->{haplotype});
    for ( my $i=0; $i<@hpl; $i++ ) {
        push @flip_var, $i unless $hpl[$i];
    }

    return unless @flip_var;
    
    my @vpos = $self->variant_positions;
    my @flip_pos = @vpos[@flip_var];

    #my $tmp = $self->{haplotype};
    #$tmp =~ tr/10/ */;
    #my @o = map{$_->[0]}@{$self->{vars}}[@vpos];
    #my @d = split(//, $tmp);
    
    foreach ( @flip_pos ) {
        @{$self->{vars}[$_]}[0,1] = @{$self->{vars}[$_]}[1,0];
        @{$self->{freqs}[$_]}[0,1] = @{$self->{vars}[$_]}[1,0];
        @{$self->{probs}[$_]}[0,1] = @{$self->{probs}[$_]}[1,0];
    }

    #my @f = map{$_->[0]}@{$self->{vars}}[@vpos];
    #for (my $i=0; $i<@f; $i++) {
    #    print $o[$i], $d[$i], $f[$i],"\n";
    #}

}


##----------------------------------------------------------------------------##

sub Sam::Seq::compute_local_haplotype{
    my $self = shift;
    my %p = (@_);
    my $frags = $p{frags};

    #print "\ncompute_local_haplotype(seed_idx => $p{seed_idx})\n";

    # $L->debug(stringify_frags(frags => $frags, max_len => length($v1), use_hpl => 1, label => "[local frags]"));

    # seed
    my $hpl = {%{$frags->[$p{seed_idx}]}};
    $hpl->{hplc} = $hpl->{hpl};
    $hpl->{hplc} =~ tr/10/01/;
    $hpl->{1} = [];
    $hpl->{0} = [];
    my $f_max_len = length($hpl->{seq});

    # use max ref similar (most "1") complement as h[0];
    if ($hpl->{hpl} =~ tr/1// < $hpl->{hplc} =~ tr/1//){
        ($hpl->{hpl}, $hpl->{hplc}) = ($hpl->{hplc}, $hpl->{hpl});
    }

    if (@$frags == 1) { # short-circuit singletons
        return $hpl;
    }

    #stringify_frags(frags => [$hpl], max_len => $self->len, use_hpl => 1, label => "[seed frag]");

    extend_suf(hpl => $hpl, frags => $frags, last_idx => $p{seed_idx});
    ## $L->debug(stringify_frags(frags => [$hpl], max_len => length($v1), use_hpl => 1); #, label => "[sufx frag]"));

    my $max_look_behind = length($frags->[$p{seed_idx}]{seq}) - Sam::Phase->FragMinOvl;
    extend_pre(hpl => $hpl, frags => $frags, first_idx => $p{seed_idx}, max_look_behind => $max_look_behind);
    ## $L->debug(stringify_frags(frags => [$hpl], max_len => length($v1), use_hpl => 1); #, label => "[prex frag]"));

    return $hpl;
}



sub extend_pre{
    my %p = (reuse_frags => [], @_);
    my $hpl = $p{hpl};
    my $f = $p{frags};

    # look back x frags for best ovl
    my $i;
    my $max_ovl = 0;
    my $min_ovl_i = $p{first_idx};
    for ($i=$p{first_idx}-1; $i>=0; $i--) {
        # can't use ovl in look-behind
        # --0000000------ < valid but missed
        # --1111---------
        # --11111-------- < no ovl
        # -------1111111-
        # my $ovl =  $f->[$i]->{end} - $hpl->{pos} + 1; # min ovl
        # last if $ovl < Sam::Phase->FragMinOvl;

        my $ovl = $hpl->{end} - $f->[$i]->{pos} + 1; # min ovl
        unless ( $ovl < Sam::Phase->FragMinOvl ){ # seen a frag with valid ovl
            $min_ovl_i = $i;
        }

        last if $f->[$i]{pos} < $hpl->{pos} - $p{max_look_behind};
    }
    $i = $min_ovl_i;

    my $x_all = $p{reuse_frags};
    if ( $i < $p{first_idx}) { # new extender
        unshift @$x_all, @{$f}[$i .. $p{first_idx}-1];
    }
    return unless @$x_all;

    # assign scores to all frags
    score_and_flip_frags(frags => $x_all, hpl => $hpl);

    # potential extenders
    my $x_pot = [];
    my $x_pre = []; # look-behind hits with too short overlap

    # remove containees, record support/coverage
    foreach ( @$x_all ) {
        if ( $_->{pos} < $hpl->{pos} ) {
            if ( $_->{end} - $hpl->{pos} + 1 >= Sam::Phase->FragMinOvl ){
                push @$x_pot, $_;
            }else {
                push @$x_pre, $_;
            }
        } else {
            add_frag_to_support_matrix($_);
        }
    }
    ## $L->debug(stringify_frags(frags => $x_con, use_hpl => 1, label => "[containee frags]"));
    return unless @$x_pot;

    # best extender (highest score)
    my $x_best;
    ($x_best, @$x_pot) = @{$x_pot}[sort { $x_pot->[$b]{score} <=> $x_pot->[$a]{score} || length($x_pot->[$b]{seq}) <=> length($x_pot->[$a]{seq}) }(0..$#$x_pot)];

    push @$x_pot, @$x_pre;
    # $L->debug(stringify_frags(frags => [$x_best], use_hpl => 1, use_scores => 1, label => "[pre-xtend best]", max_len => length($v1)));
    # $L->debug(stringify_frags(frags => $x_pot , use_hpl => 1, use_scores => 1, label => "[pre-xtend sub]", max_len => length($v1)));

    add_frag_to_support_matrix($x_best);

    my $x_len = $hpl->{pos} - $x_best->{pos};  # extension length
    my $x_hpl = substr($x_best->{hpl}, 0, $x_len);
    my $x_hplc = $x_hpl;
    $x_hplc =~ tr/10/01/;

    $hpl->{pos} = $x_best->{pos}; # set new pos
    $hpl->{hpl} = $x_hpl.$hpl->{hpl};
    $hpl->{hplc} = $x_hplc.$hpl->{hplc};


    ## $L->debug(stringify_frags(frags => [$hpl], use_hpl => 1, use_scores => 1, label => "[extended hpl]", max_len => length($v1)));

    extend_pre(hpl => $hpl, frags => $p{frags}, first_idx => $i, reuse_frags => $x_pot, max_look_behind => $p{max_look_behind});
}

sub extend_suf{
    my %p = (reuse_frags => [], @_);
    my $hpl = $p{hpl};
    my $f = $p{frags};

    # look at the next x frags for best ovl
    my $i;
    my $max_ovl = 0;
    for ($i=$p{last_idx}+1; $i<@$f; $i++) {
        # filter prior:
        # next if length($f[$i]->{seq}) < FRAG_MIN_LEN; # ignore shorties

        my $ovl = $hpl->{end} - $f->[$i]{pos} + 1; # min ovl
        last if $ovl < Sam::Phase->FragMinOvl;

        # probably insignificant pros causes trouble with containees
        #  111111
        #  11111      # strong containee
        #     111111  # week extendee
        #if ($ovl >= $max_ovl ){
        #    $max_ovl = $ovl
        #}else {
        #    last if $ovl < FRAG_MIN_MAX_OVL_FRAC * $max_ovl;
        #}
    }
    $i--;

    my $x_all = $p{reuse_frags};
    if ( $i > $p{last_idx}) { # new extenders
        push @$x_all, @{$f}[$p{last_idx}+1 .. $i ];
    }
    return unless @$x_all;

    score_and_flip_frags(frags => $x_all, hpl => $hpl);

    # potential extenders
    my $x_pot = [];
    my $x_con = [];

    # remove containees
    foreach ( @$x_all ) {
        if ( $_->{end} > $hpl->{end} ){
            push @$x_pot, $_;
        }else {
            add_frag_to_support_matrix($_);
        }
    }

    ## $L->debug(stringify_frags(frags => $x_con, use_hpl => 1, label => "[containee frags]"));
    return unless @$x_pot;

    my $x_best;
    ($x_best, @$x_pot) = @{$x_pot}[sort { $x_pot->[$b]{score} <=> $x_pot->[$a]{score} || length($x_pot->[$b]{seq}) <=> length($x_pot->[$a]{seq}) }(0..$#$x_pot)];

    add_frag_to_support_matrix($x_best);

    ## $L->debug(stringify_frags(frags => [$x_best], use_hpl => 1, use_scores => 1, label => "[suf-xtend best]", max_len => length($v1)));
    ## $L->debug(stringify_frags(frags => $x_pot, use_hpl => 1, use_scores => 1, label => "[suf-xtend sub]", max_len => length($v1)));

    my $x_len = $x_best->{end} - $hpl->{end}; # extension length
    my $x_hpl = substr($x_best->{hpl}, -$x_len);
    my $x_hplc = $x_hpl;
    $x_hplc =~ tr/10/01/;

    $hpl->{end} = $x_best->{end}; # set new end
    $hpl->{hpl}.= $x_hpl;
    $hpl->{hplc}.= $x_hplc;

    ## $L->debug(stringify_frags(frags => [$hpl], use_hpl => 1, use_scores => 1, label => "[extended hpl]", max_len => length($v1))=;

    extend_suf(hpl => $hpl, frags => $p{frags}, last_idx => $i, reuse_frags => $x_pot);
}


sub score_and_flip_frags{
    my %p = (@_);
    foreach (@{$p{frags}}) {
        my $shift = $p{hpl}{pos} - $_->{pos};
        my ($m1, $m2);
        if ($shift > 0){ # pre
            my $ovl = $_->{end} - $p{hpl}{pos}+1;
            #print "prefix\n",substr($_->{hpl}, -$ovl),"\n",$p{hpl}{hpl},"\n";
            my $s1 = substr($_->{hpl}, -$ovl) ^ $p{hpl}{hpl};
            $m1 = $s1 =~ tr/\0//; # count matches
            my $s2 = substr($_->{hpl}, -$ovl) ^ $p{hpl}{hplc};
            $m2 = $s2 =~ tr/\0//; # count matches
        } else { # suf or containee
            #print "suffix\n",substr($p{hpl}{hpl}, abs($shift)),"\n",$_->{hpl},"\n";
            my $s1 = substr($p{hpl}{hpl}, abs($shift)) ^ $_->{hpl};
            $m1 = $s1 =~ tr/\0//; # count matches
            my $s2 = substr($p{hpl}{hplc}, abs($shift)) ^ $_->{hpl};
            $m2 = $s2 =~ tr/\0//; # count matches
        }

        if($m2 < $m1){
            $_->{score} =  $m1;
            $_->{flipped} = 0;
        }else {
            $_->{score} = $m2;
            $_->{hpl} =~ tr/10/01/;
            $_->{flipped} = 1;
        }
    }
}

sub add_frag_to_support_matrix{
    my ($f) = (@_);
    my $k = $f->{pos};
    foreach my $state ( split(//, $f->{hpl})) {
        $support_matrix->{$state}[$k]++;
        $k++;
    }
}

sub print_haplotypes{
    my %p = (@_);
    my $suf = length($p{ref}) - ($p{pos} + length($p{hpl}[0]));
    print "HAPLOTYPES:\n";
    print "-" x $p{pos} , $p{hpl}[0], "-" x $suf,"\n";
    print "-" x $p{pos} , $p{hpl}[1], "-" x $suf,"\n";
};

sub stringify_frags{
    #return unless $debug;
    my %p = (@_);
    my @re;
    push @re, $p{label},"\n" if $p{label};
    foreach (@{$p{frags}}) {
        my $k = $p{use_hpl} ? "hpl" : "seq";
        my $suf = $p{max_len} ? $p{max_len} - ($_->{pos} + length($_->{$k})) : 0;
        my $score = "";
        if ( $p{use_scores} && exists $_->{score} ) {
            $score = " [$_->{score}]". ($_->{flipped} ? " !" : "");
        }
        push @re, "-" x $_->{pos}. $_->{$k}. "-" x $suf . $score . "\n";
    }
    return @re;
}


##----------------------------------------------------------------------------##
## Closures

sub _init_closures{
    no strict 'refs';

    foreach my $attr ( keys %ATTR_CLASS ){
        next if __PACKAGE__->can($attr); # don't overwrite explicitly declared methods
        *{__PACKAGE__ . "::$attr"} = sub {
            die (((caller 0)[3])." expects a single value, got @_\n") if @_>2;
            $ATTR_CLASS{$attr} = $_[1] if @_ == 2;
            return $ATTR_CLASS{$attr};
        }
    }

    # might wann treat ARRAYs and HASHs differently
}


__PACKAGE__->_init_closures(); # close to eof


=head1 AUTHOR

Thomas Hackl, E<lt>thackl@lim4.deE<gt>

=cut

1;
