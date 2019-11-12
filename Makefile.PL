use ExtUtils::MakeMaker;

WriteMakefile(NAME      => 'RyanLabPhylogenomicTools',
              AUTHOR    => 'Joseph Ryan <joseph.ryan@whitney.ufl.edu',
              PREREQ_PM => { 'JFR::Fasta' => '1.04' },
              EXE_FILES => ['hmm2aln.pl','stockholm2fasta.pl']);
