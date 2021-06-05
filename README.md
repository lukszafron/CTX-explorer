# CTX-explorer

Welcome to the CTX-Explorer app. 

This open-source Python 3 program was developed for identification of inter- and intrachromosomal translocations.

Additional Python libraries that need to be installed for this app to work correcly:
subprocess, re, os, gzip, sys, getopt, statistics, itertools, termcolor

Additional software that needs to be installed for this app to work correctly:
samtools, bamtools

The following options are available:

        -h, --help:             prints this help message.
        -t, --threads:          defines the number of CPU threads that are to be used (default: 1).
        -T, --tmpdir:           defines the directory where temporary files are to be stored (default: current dir).
        -q, --qual:             defines the minimal Phred quality score of mappings to be used (default: 20).
        -b, --bamfile:          the name of a BAM file to be used (MANDATORY).
        -I, --intractx:         indicates if intra-chromosomal translocations should be evaluated (default: 'no'). This option significantly increases memory usage.
        -1, --chrom1:           the name of the first chromosome (optional) involved in a translocation (e.g., 'chr1').
        -2, --chrom2:           the name of the second chromosome (optional) involved in a translocation (e.g., 'chr2').
        -l, --tlen:             minimal distance (in nucleotides) between the first and the second read forming a read pair (default: 1,000,000).
        -i, --insert:           the maximum insert size used for identification of CTX-supporting reads and hit pairs (default: median insert size * 2).
        -n, --nohits:           minimal number of hits for a translocation to be stored in the final report (default: 2).
        -N, --nohits_sec:       minimal number of hits per the second chromosome for a translocation to be stored in the final report (default: 2).
        -s, --min_size:         minimal size of a hit group for a translocation to be stored in the final report (default: 2).
        -d, --no_filter:        specifies whether the filtering by the number of supporting reads should be turned off (default: 'no').
        -p, --prefix:           specifies a prefix of a file in which the final report will be saved (default: 'output').
        -g, --gzipped:          indicates whether the final report file should be gzipped (default: 'no').
        -v, --version:          prints the version of this program.

