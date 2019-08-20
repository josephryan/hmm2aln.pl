# hmm2aln.pl
script to automate gene family phylogeny based using HMM

This repo is under construction.

###  --nofillcnf
HMMer can make alignments from an hmmsearch but often misses the ends of alignments. hmm2aln.pl rescues these missing bits, and uses the --nofillcnf option (see nofill.hox.conf file) to make sure it doesn't rescue garbage.  The option requires that certain residue positions only be rescued if they are certain amino acids.

