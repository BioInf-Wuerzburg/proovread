package Fasta::Store;

use strict; 
use warnings;

use Log::Log4perl;

use Handy qw(require_exe);
use Fasta::Seq;

our $VERSION = '0.1.0';

=head1 NAME

Fasta::Store 

=head1 SYNOPSIS

  use Fasta::Store;

=head1 DESCRIPTION

Random region access to index fasta files using samtools faidx.

=cut

my @ATTR_SCALAR = qw(file samtools_path samtools fai fai_cmd);

my %SELF;
@SELF{@ATTR_SCALAR} = (undef) x scalar @ATTR_SCALAR;

my $L = Log::Log4perl->get_logger();

=head1 METHODS

=head2 new

=cut

sub new{
    my $class = shift;
    my $self = {
        %SELF,
        samtools => 'samtools',
        @_
    };

    bless $self, $class;

    die (((caller 0)[3]).": file required\n") unless defined ($self->file) && length($self->file);
    die (((caller 0)[3]).": ".$self->file." doesn't exist or isn't accessible\n") unless -e $self->file && -r _;

    $self->samtools(require_exe(join("/", grep{$_}($self->samtools_path, $self->samtools))));

    $self->fai($self->file.".fai");
    $self->fai_cmd($self->samtools." faidx ".$self->file);
    unless (-e $self->fai) {
        $L->debug($self->fai_cmd);
        my $fai_cmd = $self->fai_cmd;
        qx($fai_cmd);
        $? && die $?;
    }
    
    return $self;
}

=head2 file, fai, samtools, samtools_path

Get/set ...

=cut

sub _init_accessors{
    no strict 'refs';

    # generate simple accessors closure style
    foreach my $attr ( @ATTR_SCALAR ) {
        next if $_[0]->can($attr); # don't overwrite explicitly created subs
        *{__PACKAGE__ . "::$attr"} = sub {
            $_[0]->{$attr} = $_[1] if @_ == 2;
            return $_[0]->{$attr};
        }
    }
}

=head2 fetch

=cut

sub fetch{
    my ($self, $id) = @_;
    die (((caller 0)[3]).": id or region required\n") unless defined $id;
    my $cmd = $self->fai_cmd." $id";
    my $fa = scalar qx/$cmd/;
    $? && die $?>>8;
    return Fasta::Seq->new($fa);
}

__PACKAGE__->_init_accessors();

=head1 AUTHOR

Thomas Hackl, E<lt>thackl@lim4.deE<gt>

=cut

1;

