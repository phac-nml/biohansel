Command-Line
============

If you run ``hansel -h``, you should see the following usage statement:

.. code-block::

    usage: hansel [-h] [-s SCHEME] [--scheme-name SCHEME_NAME]
                  [-p forward_reads reverse_reads] [-i fasta_path genome_name]
                  [-D INPUT_DIRECTORY] [-o OUTPUT_SUMMARY]
                  [-O OUTPUT_TILE_RESULTS] [-S OUTPUT_SIMPLE_SUMMARY] [--force]
                  [--json] [--min-kmer-freq MIN_KMER_FREQ]
                  [--max-kmer-freq MAX_KMER_FREQ]
                  [--low-cov-depth-freq LOW_COV_DEPTH_FREQ]
                  [--max-missing-tiles MAX_MISSING_TILES]
                  [--min-ambiguous-tiles MIN_AMBIGUOUS_TILES]
                  [--low-cov-warning LOW_COV_WARNING]
                  [--max-intermediate-tiles MAX_INTERMEDIATE_TILES] [-t THREADS]
                  [-v] [-V]
                  [F [F ...]]

    Subtype microbial genomes using SNV targeting k-mer subtyping schemes.
    Includes schemes for Salmonella enterica spp. enterica serovar Heidelberg and Enteritidis subtyping.
    Developed by Geneviève Labbé, James Robertson, Peter Kruczkiewicz, Marisa Rankin, Matthew Gopez, Chad R. Laing, Philip Mabon, Kim Ziebell, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Philip Mabon, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Marisa Rankin, Joanne MacKinnon, Roger Johnson, John H.E. Nash.

    positional arguments:
      F                     Input genome FASTA/FASTQ files (can be Gzipped)

    optional arguments:
      -h, --help            show this help message and exit
      -s SCHEME, --scheme SCHEME
                            Scheme to use for subtyping (built-in: "heidelberg",
                            "enteritidis"; OR user-specified:
                            /path/to/user/scheme)
      --scheme-name SCHEME_NAME
                            Custom user-specified SNP substyping scheme name
      -p forward_reads reverse_reads, --paired-reads forward_reads reverse_reads
                            FASTQ paired-end reads
      -i fasta_path genome_name, --input-fasta-genome-name fasta_path genome_name
                            fasta file path to genome name pair
      -D INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                            directory of input fasta files (.fasta|.fa|.fna) or
                            FASTQ files (paired FASTQ should have same basename
                            with "_\d\.(fastq|fq)" postfix to be automatically
                            paired) (files can be Gzipped)
      -o OUTPUT_SUMMARY, --output-summary OUTPUT_SUMMARY
                            Subtyping summary output path (tab-delimited)
      -O OUTPUT_TILE_RESULTS, --output-tile-results OUTPUT_TILE_RESULTS
                            Subtyping tile matching output path (tab-delimited)
      -S OUTPUT_SIMPLE_SUMMARY, --output-simple-summary OUTPUT_SIMPLE_SUMMARY
                            Subtyping simple summary output path
      --force               Force existing output files to be overwritten
      --json                Output JSON representation of output files
      --min-kmer-freq MIN_KMER_FREQ
                            Min k-mer freq/coverage
      --max-kmer-freq MAX_KMER_FREQ
                            Max k-mer freq/coverage
      --low-cov-depth-freq LOW_COV_DEPTH_FREQ
                            Frequencies below this coverage are considered low
                            coverage
      --max-missing-tiles MAX_MISSING_TILES
                            Decimal proportion of maximum allowable missing tiles
                            before being considered an error. (0.0 - 1.0)
      --min-ambiguous-tiles MIN_AMBIGUOUS_TILES
                            Minimum number of missing tiles to be considered an
                            ambiguous result
      --low-cov-warning LOW_COV_WARNING
                            Overall tile coverage below this value will trigger a
                            low coverage warning
      --max-intermediate-tiles MAX_INTERMEDIATE_TILES
                            Decimal proportion of maximum allowable missing tiles
                            to be considered an intermediate subtype. (0.0 - 1.0)
      -t THREADS, --threads THREADS
                            Number of parallel threads to run analysis (default=1)
      -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
                            show debug info)
      -V, --version         show program's version number and exit
