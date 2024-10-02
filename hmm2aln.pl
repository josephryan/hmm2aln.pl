#!perl

$|++;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use JFR::Fasta;
use Data::Dumper;

our $VERSION = '2.0.0';

# requires hmmsearch from hmmer package http://hmmer.org/
our $HMMSEARCH = 'hmmsearch';
our $ESL_REFORMAT = 'esl-reformat';

MAIN: {
    my $rh_o = get_opts();
    my $ra_fa   = get_fasta($rh_o);
    $rh_o->{'count'} = 0;
    my $ra_nfc = get_nfc($rh_o->{'nofillcnf'}) if ($rh_o->{'nofillcnf'});
    FA: foreach my $fa (@{$ra_fa}) {
        my ($rfa,$ra_orig_defs) = rename_defs($rh_o,$fa); #00 files
        my $pre = $rh_o->{'name'};
        run_hmmsearch($rfa,$rh_o);  #01x files
        if (-z "$pre/01b-$rh_o->{'count'}-$pre.aln") {
            warn "domain not found in $fa\n";
            next FA;
        }
        stockholm2fasta($rh_o);     #02x files
        remove_gaps($rh_o);         #03 file

        unless (-z $rfa) {  # unless file is empty
            my $fa_aln = deal_with_gaps($rh_o,$rfa,$ra_orig_defs); 
            if ($rh_o->{'nofillcnf'}) {
                replace_w_unfilled($rh_o,\$fa_aln,$ra_nfc,$ra_orig_defs);
            }
            my $final = number_multi_domains($fa_aln);
            print $final;
        }
        $rh_o->{'count'}++;
    }
    clean_up($rh_o) unless ($rh_o->{'no_clean'});
}

sub replace_w_unfilled {
    my $rh_o = shift;
    my $rs_fa_aln = shift;
    my $ra_nfc = shift;
    my $ra_defs = shift;
    my %remove_fill = ();

    my @lines = split /\n/, ${$rs_fa_aln};
    foreach my $ra_crit (@{$ra_nfc}) {
        get_nofill_ids(\@lines,$ra_crit,\%remove_fill);
    }
    my $rh_unfilled = get_unfilled_seqs($rh_o,$ra_defs);
    ${$rs_fa_aln} = '';
    for (my $i = 0; $i < @lines; $i += 2) {
        my $def = $lines[$i];
        my $seq = $lines[$i+1];
        ${$rs_fa_aln} .= $def . "\n";
        if ($remove_fill{$def}) {
            ${$rs_fa_aln} .= $rh_unfilled->{$def} . "\n";
        } else {
            ${$rs_fa_aln} .= $seq . "\n";
        }
    }
}

sub get_unfilled_seqs {
    my $rh_o = shift;
    my $ra_defs = shift;
    my %seqs = ();
    my $nogaps_fa = "$rh_o->{'name'}/03-$rh_o->{'count'}-$rh_o->{'name'}.fa.nogaps";
    my $fp = JFR::Fasta->new($nogaps_fa);
    while (my $rec = $fp->get_record()) {
        $rec->{'def'} =~ m/^>(\d+)\/\d+-\d+/ or die "unexpected";
        my $def = $ra_defs->[$1];
        $seqs{$def} = $rec->{'seq'};
    }
    return \%seqs;
}

sub get_nofill_ids {
    my $ra_lines = shift;
    my $ra_c = shift;
    my $rh_rf = shift;

#    my $fp = JFR::Fasta->new($file);
    for (my $i = 0; $i < @{$ra_lines}; $i += 2) {
        my $missed = 0;
        my $def = $ra_lines->[$i];
        my $seq = $ra_lines->[$i+1];
        foreach my $ra_data (@{$ra_c}) {
            my $sub = substr($seq,$ra_data->[1],1);
            my $char_class = '[' . $ra_data->[2] . ']';
            $missed++ unless ($sub =~ m/^$char_class$/i);
        }
        $rh_rf->{$def}++ if ($missed);
    }
}

sub get_no_f {
    my $file = shift;
    my $rh_rf = shift;
    my $fp = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        my $nsub = substr($rec->{'seq'},19,1);
        next if ($nsub =~ m/[-FYH]/);
        $rh_rf->{$rec->{'def'}}++;
    }
}

# routine to parse and check formatting of --nofillcnf
sub get_nfc {
    my $conf = shift;
    my @nfc = ();
    open IN, $conf or die "cannot open $conf:$!";
    while (my $line = <IN>) {
        next if ($line =~ m/^\s*$/);
        next if ($line =~ m/^\s*#/);
        $line =~ s/\s+$//;
        my @f = split /\s+/, $line;
        die "error in format of $conf" unless (scalar(@f) == 3);
        die "error in format of $conf" unless ($f[0] =~ m/^\d+$/);
        die "error in format of $conf" unless ($f[1] =~ m/^\d+$/);
        $f[1]--;
        $f[0]--;
        die "unexpected val in $conf" unless (($f[1] >= 0) && ($f[0] >= 0));
        push @{$nfc[$f[0]]}, \@f;
    }
    return \@nfc;
}

sub number_multi_domains {
    my $fa_str = shift;
    my @lines = split /\n/, $fa_str;
    my %defs = ();
    my %counter = ();
    my $new_str = '';
    foreach my $l (@lines) {
        $defs{$l}++ if ($l =~ m/^>/);
    }
    foreach my $l (@lines) {
        if ($defs{$l} && $defs{$l} > 1) {
            $counter{$l}++;
            $new_str .= "$l.$counter{$l}\n";
        } else {
            $new_str .= "$l\n";
        }
    }
    return $new_str;
}

sub clean_up {
    my $rh_o = shift;
    my $n = $rh_o->{'name'};
    my @files = ("00-$n.renamed.fa", "01a-hmmsearch.out", "01b-$n.aln",
                 "01c-$n.hmmsearch", "02a-esl-reformat.out", "02b-$n.fa",
                 "03-$n.fa.nogaps");
    my $err = 0;
    foreach my $f (@files) {
        unlink "$n/$f" || $err++;
    }
    rmdir $n unless ($err);
    warn "some intermediate files in $n/ were not deleted\n" if ($err);
}

sub deal_with_gaps {
    my $rh_o    = shift;
    my $fa      = shift;
    my $ra_defs = shift;
    my $ra_seqs = _get_seqs($fa);
    my $fa_aln  = '';

    my $nogaps_fa = "$rh_o->{'name'}/03-$rh_o->{'count'}-$rh_o->{'name'}.fa.nogaps";

    my $fp = JFR::Fasta->new($nogaps_fa);

    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        $id =~ m/^(\d+)\/(\d+)-(\d+)/;
        my $i = $1; 
        my $start = $2;
        my $end = $3;
        if ($rh_o->{'fillendgaps'} && 
           ($rec->{'seq'} =~ m/^-/ || $rec->{'seq'} =~ m/-$/) ) {
            if ($rec->{'seq'} =~ m/^(-+)/) {
                my $s_gaps = length $1;
                my $j = $start - $s_gaps - 1;
                if ($j >= 0) {
                    my $replace = substr $ra_seqs->[$i]->{'seq'},$j,$s_gaps;
                    $rec->{'seq'} =~ s/^-{$s_gaps}/$replace/;
                }
            }
            if ($rec->{'seq'} =~ m/(-+)$/) {
                my $e_gaps = length $1;
                my $len = length $ra_seqs->[$i]->{'seq'};
                if ($len - $end +1 - $e_gaps > 0) {
                    my $replace = substr $ra_seqs->[$i]->{'seq'},$end,$e_gaps;
                    $rec->{'seq'} =~ s/-{$e_gaps}$/$replace/;
                }
            }
            
        } 
        $fa_aln .= "$ra_defs->[$i]\n$rec->{'seq'}\n";
    }
    return $fa_aln;
}
    

sub _get_seqs {
     my $file = shift;
     my @seqs = ();
     my $fp = JFR::Fasta->new($file);
     while (my $rec = $fp->get_record()) {
         $rec->{'seq'} =~ s/\.//g;
         push @seqs, $rec;
     }
     return \@seqs;
}


sub remove_gaps {
    my $rh_o = shift;
    my %seqs = ();
    my $ra_gaps = _get_gap_positions($rh_o,\%seqs);
    _print_sequences($rh_o,$ra_gaps,\%seqs);
}

sub _print_sequences {
    my $rh_o = shift;
    my $ra_gaps = shift;
    my $rh_seqs = shift;
    my $out = "$rh_o->{'name'}/03-$rh_o->{'count'}-$rh_o->{'name'}.fa.nogaps";
    open OUT, ">$out" or die "cannot open >$out:$!";
    foreach my $def (keys %{$rh_seqs}) {
        print OUT "$def\n";
        for (my $i=0; $i < @{$rh_seqs->{$def}}; $i++) {
            print OUT "$rh_seqs->{$def}->[$i]" unless ($ra_gaps->[$i]);
        }
        print OUT "\n";
    }
    close OUT;
}

sub _get_gap_positions {
    my $rh_o = shift;
    my $rh_seqs = shift;
    my @gaps = ();
    my $fp = JFR::Fasta->new("$rh_o->{'name'}/02b-$rh_o->{'count'}-$rh_o->{'name'}.fa");
    while (my $rec = $fp->get_record()) {
        my @aas = split /|/, $rec->{'seq'};
        $rh_seqs->{$rec->{'def'}} = \@aas;
        for (my $i = 0; $i < @aas; $i++) {
            $gaps[$i]++ if ($aas[$i] =~ m/[a-z*]/);
        }
    }
    return \@gaps;
}

sub stockholm2fasta {
    my $rh_o = shift;
    my $pre = $rh_o->{'name'};
    my $cmd = "$ESL_REFORMAT -o $pre/02b-$rh_o->{'count'}-$pre.fa afa ";
    $cmd .= "$pre/01b-$rh_o->{'count'}-$pre.aln ";
    $cmd .= "> $pre/02a-$rh_o->{'count'}-esl-reformat.out 2>&1";
    system $cmd;
}

sub run_hmmsearch {
    my $fa = shift;
    my $rh_o = shift;
    
    my $dir = $rh_o->{'name'};
    my $cmd = "hmmsearch ";
    $cmd .= "--cpu $rh_o->{'threads'} " if ($rh_o->{'threads'});
    $cmd .= "-A $dir/01b-$rh_o->{'count'}-$rh_o->{'name'}.aln ";
    $cmd .= "-o $dir/01c-$rh_o->{'count'}-$rh_o->{'name'}.hmmsearch ";
    $cmd .= "$rh_o->{'hmm'} $fa  > $dir/01a-$rh_o->{'count'}-hmmsearch.out 2>&1";
    system $cmd;
}

sub get_fasta {
    my $rh_o = shift;
    my @fasta = ();
    if ($rh_o->{'fasta'}) {
        $fasta[0] = $rh_o->{'fasta'};
    } elsif ($rh_o->{'fasta_dir'}) {
        my $dir = $rh_o->{'fasta_dir'};
        die "$dir is not a directory" unless (-d $dir);
        opendir DIR, $dir or die "cannot opendir $dir:$!";
        my @files = grep { !/^\./ } readdir DIR;
        foreach my $f (@files) {
            if (_is_fasta("$dir/$f")) {
                push @fasta, "$dir/$f";
            }
        }
        die "no fasta files in $dir\n" unless (scalar(@fasta));
    }
    return \@fasta;
}

# tests if the first non-blank line starts with >
sub _is_fasta {
    my $file = shift;
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next if ($line =~ m/^\s*$/);
        if ($line =~ m/^\s*>/) {
            return 1;
        } else {
            return 0;
        }
    }
}

sub rename_defs {
    my $rh_o = shift;
    my $fa = shift;
    my $fp = JFR::Fasta->new($fa);
    my @orig_defs = ();
    my $count = 0;
    my $rfa = "$rh_o->{'name'}/00-$rh_o->{'name'}.renamed.fa";
    open OUT, ">$rfa" or die "cannot open >$rfa:$!";
    while (my $rec = $fp->get_record()) {
        push @orig_defs, $rec->{'def'};
        print OUT ">$count\n$rec->{'seq'}\n";
        $count++;
    }
    return ($rfa,\@orig_defs);
}


sub get_opts {
    my %opts = ('hmm' => '', 'name' => '', 'fasta' => '', 'fasta_dir' => '',
                'no_clean' => '', 'nofillcnf' => '', 'help' => '', 
                'version' => '');
    my $opt_results = Getopt::Long::GetOptions(
                                      'hmm=s' => \$opts{'hmm'},
                                      'name=s' => \$opts{'name'},
                                      'fasta=s' => \$opts{'fasta'},
                                      'fasta_dir=s' => \$opts{'fasta_dir'},
                                      'threads=i' => \$opts{'threads'},
                                      'no_clean' => \$opts{'no_clean'},
                                      'fillendgaps' => \$opts{'fillendgaps'},
                                      'nofillcnf=s' => \$opts{'nofillcnf'},
                                      'help' => \$opts{'help'},
                                      'version' => \$opts{'version'});
    die "$0 version $VERSION\n" if ($opts{'version'});
    pod2usage({-exitval => 0, -verbose => 2}) if($opts{'help'});
    usage() unless ($opts{'hmm'} && $opts{'name'});
    usage() unless ($opts{'fasta'} || $opts{'fasta_dir'});
    if ($opts{'fasta'} && $opts{'fasta_dir'}) {
        die "can not provide both --fasta AND --fasta_dir\n";
    }
    die "error: directory \"$opts{'name'}\" exists\n" if (-d $opts{'name'});
    die "must use --fillendgaps when using --nofillcnf" if (!$opts{'fillendgaps'} && $opts{'nofillcnf'});
    mkdir $opts{'name'} or die "cannot mkdir $opts{'name'}:$!";
    return \%opts;    
}

sub usage {
    die "usage: $0 --hmm=<hmmfile> --name=<name> {--fasta=<fasta>|--fasta_dir=<fasta_dir>} [--threads=<num>] [--no_clean] [--fillendgaps] [--nofillcnf=<nofill.conf>] [--help] [--version]\n";
}

__END__

=head1 NAME

B<hmm2align.pl> - use hmm to build alignment

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

hmm2aln.pl --hmm=<hmmfile> --name=<name> {--fasta=<fasta>|--fasta_dir=<fasta_dir>} [--threads=<num>] [--no_clean] [--fillendgaps] [--nofillcnf=<nofill.conf>] [--help] [--version]

=head1 OPTIONS

=over

=item B<--hmm>

hmmfile (e.g. download http://pfam.xfam.org/family/PF00096/hmm)

=item B<--name>
name used for output files

=item B<--fasta>
FASTA-formatted file with multiple sequences

=item B<--fasta_dir>
directory with FASTA-formatted files

=item B<--threads>
number of threads to use for search

=item B<--no_clean>
do not remove intermediate files created by hmmer

=item B<--fillendgaps>
In versions 1.0.1 and earlier, this was the default behavior. When making alignments from an hmmsearch, HMMer often misses the ends of alignments. When this parameter is added to the command line, hmm2aln.pl will pad end gap positions with those amino acids adjacent to the match. For example, if a domain was detected starting at position 103 in an amino acid sequence, and this position matched position 3 of the HMM, positions 101 and 102 from the sequence would be included in the recovered domain (instead of 2 gaps in those positions). This can be controled with the --nofillcnf option, which requires that fills only occur when certain residues are present at certain positions.

=item B<--nofillcnf>
HMMer can make alignments from an hmmsearch but often misses the ends of alignments. hmm2aln.pl rescues these missing bits, and uses the --nofillcnf option (see nofill.hox.conf file) to make sure it doesn't rescue garbage.  The option requires that certain residue positions only be rescued if they are certain amino acids. See example file (https://github.com/josephryan/hmm2aln.pl/blob/master/nofill.hox.conf) for format.

=item B<--help>
print this manual

=item B<--version>
print the version. Overrides all other options.

=head1 DESCRIPTION

This program takes an HMM and a multi-sequence FASTA file (or directory of FASTA files) and produces a multi-sequence alignment.

=head1 COPYRIGHT

Copyright (C) 2018,2019,2020,2021,2022,2023,2024 Joseph F. Ryan

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
