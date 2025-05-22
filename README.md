# hmm2aln.pl
script to automate gene family phylogeny based using HMM

## NOTE: 2.0 version (and beyond) does not fill end gaps by default. To fill end gaps use --filendgaps option

### REQUIREMENTS

* HMMER
http://hmmer.org/download.html

* JFR::Fasta
https://github.com/josephryan/JFR-PerlModules

### INSTALL

```
   git clone https://github.com/josephryan/hmm2aln.pl
   cd hmm2aln.pl
```

To install these modules and scripts type the following:

```
   perl Makefile.PL
   make
   sudo make install
```
   
To install without root privelages try:

```
   perl Makefile.PL PREFIX=/home/myuser/scripts
   make
   make install
```

Or just copy scripts to a directory and run them like this:

```
   perl hmm2aln.pl
```

## hmm2aln.pl 

script takes an HMM and a multi-sequence FASTA file (or directory of FASTA files) and produces an alignment.

### SYNOPSIS

hmm2aln.pl --hmm=<hmmfile> --name=<name> {--fasta=<fasta>|--fasta_dir=<fasta_dir>} [--threads=<num>] [--no_clean] [--fillendgaps] [--nofillcnf=<nofill.conf>] [--help] [--version]

### OPTIONS

       --hmm
            hmmfile (e.g. download http://pfam.xfam.org/family/PF00096/hmm)

       --name
            name used for output files

       --fasta
            FASTA-formatted file with multiple sequences

       --fasta_dir
            directory with FASTA-formatted files

       --threads
            number of threads to use for search

       --no_clean
            do not remove intermediate files created by hmmer

    --fillendgaps In versions 1.0.1 and earlier, this was the default behavior. When making alignments from an hmmsearch, HMMer often misses the ends of alignments. When this parameter is added to the command line, hmm2aln.pl will pad end gap positions with those amino acids adjacent to the match. For example, if a domain was detected starting at position 103 in an amino acid sequence, and this position matched position 3 of the HMM, positions 101 and 102 from the sequence would be included in the recovered domain (instead of 2 gaps in those positions). This can be controled with the --nofillcnf option, which requires that fills only occur when certain residues are present at certain positions.

    --nofillcnf HMMer can make alignments from an hmmsearch but often misses the ends of alignments. hmm2aln.pl rescues these missing bits, and uses the --nofillcnf option (see nofill.hox.conf file) to make sure it doesn't rescue garbage. The option requires that certain residue positions only be rescued if they are certain amino acids. See example file (https://github.com/josephryan/hmm2aln.pl/blob/master/nofill.hox.conf) for format.

       --help
            print this manual

       --version
            print the version. Overrides all other options.

## stockholm2fasta.pl

script used by hmm2aln.pl to convert Stockholm-formatted multiple sequence alignment to FASTA-formatted sequence alignment. Uses esl-reformat program, which comes with HMMer, but since esl-reformat will not print fasta with gaps this script uses esl-reformat to print a fasta (long deflines) and phys format (short deflines but seqs w/gaps) and then prints a hybrid of the 2 (long defs and seqs w/gaps)

## Cited

hmm2aln.pl has been used in the multiple studies:
* https://doi.org/10.1093/gbe/evac144
* https://doi.org/10.1093/gbe/evac172
* https://doi.org/10.1093/molbev/msad137
* https://doi.org/10.1093/g3journal/jkaf110
* https://doi.org/10.1093/gbe/evaa060
* https://doi.org/10.1093/nargab/lqae072

## BUGS
    Please report them to Joseph Ryan <joseph.ryan@whitney.ufl.edu>

## COPYRIGHT
    Copyright (C) 2018,2019 Joseph F. Ryan

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
