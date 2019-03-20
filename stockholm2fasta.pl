#!/usr/bin/perl

# Sat Aug 25 16:14:08 EDT 2018
# 
# esl-reformat will not print fasta with gaps
# this script uses esl-reformat to print a fasta (long deflines)
# and phys format (short deflines but seqs w/gaps) and then prints
# a hybrid of the 2 (long defs and seqs w/gaps)

use strict;
use warnings; 
use JFR::Fasta;
use File::Temp qw/tmpnam/;
use Data::Dumper;

our $VERSION = 0.01;

MAIN: {
    my $fa_out = tmpnam();
    my $phy_out = tmpnam();
    my $stk = $ARGV[0] or die "usage: $0 STOCKHOLM_FORMATTED_FILE\n";
    system "esl-reformat fasta $stk > $fa_out";
    system "esl-reformat phylips $stk > $phy_out";
    my $ra_fa  = get_fa("$fa_out");
    my $ra_phy = get_phy("$phy_out");
    print_fa($ra_fa,$ra_phy);
}

sub print_fa {
    my $ra_fa = shift;
    my $ra_phy = shift;
    for (my $i = 0; $i < @{$ra_phy}; $i++) {
        print "$ra_fa->[$i]\n$ra_phy->[$i]\n";
    }
}

sub get_phy {
    my $phy  = shift;
    my @seqs = ();
    open IN, $phy or die "cannot open $phy:$!";
    my $first = <IN>;
    my $last = '';
    while (my $line = <IN>) {
        $line =~ s/\s+$//;
        my @f = split /\s+/, $line;
        if (scalar(@f) > 1) {
            push @seqs, $last if ($last);
            $last = $f[-1];
        } else {
            $last .= $f[0];
        }
    }
    push @seqs, $last;
    return \@seqs;
}
sub get_fa {
    my $fa   = shift;
    my @defs = ();
    my $fp   = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        push @defs, $rec->{'def'};
    }
    return \@defs;
}
