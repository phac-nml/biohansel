Command-Line
============

Here, the arguments needed to run biohansel effectively are displayed. The required and additional arguments are shown below to see what must be included in a run. 

Required
########

Make sure to be in the directory containing all of the data needed to run a command or that the path to the input data is put into the command following the argument.

- Subtyping Scheme

    - use -s "scheme"

- Output/Results Files (any combination so long as there is at least one specified. Details in `Output <output.html>`_)

    - use -S "filename.tab" | for tech_results.tab
    - use -o "filename.tab" | for results.tab
    - use -O "filename.tab" | for match_results.tab
    
- Input data

    - use -i <path/to/fasta> | to specify fasta file to analyze
    - use -p <path/to/forward_reads> <path/to/reverse_reads> | to analyze paired reads
    - use -D <path/to/directory> | to analyze a full directory of data into 1 file


Additional
##########

**If any of these arguments are left off of the command used to run biohansel, they will be set to default values for the given analysis.**

| 
| -M "metadata_scheme.tsv"  -->  Used to input a metadata scheme that follows all requirements
|        found in `input <input.html>`_
|
| \\--force  -->  Forces the existing outputs to be overwritten
|
| \\--json  -->  Output JSON representation of output files
|
| \\--min-kmer-freq <#>  -->  Minimum k-mer coverage needed for a raw reads fastq file to be
|       considered acceptable by the quality control module (default is 8) (cannot lower past 8 in
|       the current build)
|
| \\--max-kmer-freq <#>  -->  Maximum k-mer coverage for a raw reads fastq file to be considered
|       acceptable (default is 1000)
|
| \\--low-cov-depth-freq <#>  -->  Coverage frequencies of raw read fastq files below this value are
|       considered as low coverage (default is 20)
|
| \\--max-missing-kmers <#>  -->  Decimal proportion of maximum allowable missing kmers before being
|       considered an error (0.0 - 1.0) (default is 0.05 or 5%)
|
| \\--min-ambiguous-kmers <#>  -->  Minimum number of missing kmers to be considered an ambiguous
|       result (default is 3)
|
| \\--low-cov-warning <#>  -->  Overall kmer coverage below this value will trigger a low coverage 
|       warning on raw read fastq files. (default is 20) 
|
| \\--max-intermediate-kmers <#>  -->  Decimal proportion of maximum allowable missing kmers
|       (0.0 - 1.0) to be considered an intermediate subtype (default is 0.05)
| 
| --threads <#_CPUs>  -->  Number of parallel threads used to run the analysis (default = 1)
|
| -v  -->  Verbose: Logs verbosity levels where -v == show warnings and -vv == show debug info
|
| -V  -->  Displays the version of biohansel installed
|

Hansel Help Command
###################


If you run ``hansel -h``, you will be provided with additional information for most of the commands along with following usage statement:

.. code-block:: bash

    usage: hansel [-h] [-s SCHEME] [--scheme-name SCHEME_NAME]
                  [-p forward_reads reverse_reads] [-i fasta_path genome_name]
                  [-D INPUT_DIRECTORY] [-o OUTPUT_SUMMARY]
                  [-O OUTPUT_KMER_RESULTS] [-S OUTPUT_SIMPLE_SUMMARY] [--force]
                  [--json] [--min-kmer-freq MIN_KMER_FREQ]
                  [--max-kmer-freq MAX_KMER_FREQ]
                  [--low-cov-depth-freq LOW_COV_DEPTH_FREQ]
                  [--max-missing-kmers MAX_MISSING_KMERS]
                  [--min-ambiguous-kmers MIN_AMBIGUOUS_KMERS]
                  [--low-cov-warning LOW_COV_WARNING]
                  [--max-intermediate-kmers MAX_INTERMEDIATE_KMERS] [-t THREADS]
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
      -M SCHEME_METADATA, --scheme-metadata SCHEME_METADATA
                            Scheme subtype metadata table (CSV or tab-delimited
                            format; must contain "subtype" column)
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
      -O OUTPUT_KMER_RESULTS, --output-kmer-results OUTPUT_KMER_RESULTS
                            Subtyping kmer matching output path (tab-delimited)
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
      --max-missing-kmers MAX_MISSING_KMERS
                            Decimal proportion of maximum allowable missing kmers
                            before being considered an error. (0.0 - 1.0)
      --min-ambiguous-kmers MIN_AMBIGUOUS_KMERS
                            Minimum number of missing kmers to be considered an
                            ambiguous result
      --low-cov-warning LOW_COV_WARNING
                            Overall kmer coverage below this value will trigger a
                            low coverage warning
      --max-intermediate-kmers MAX_INTERMEDIATE_KMERS
                            Decimal proportion of maximum allowable missing kmers
                            to be considered an intermediate subtype. (0.0 - 1.0)
      -t THREADS, --threads THREADS
                            Number of parallel threads to run analysis (default=1)
      -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
                            show debug info)
      -V, --version         shows the program version number and exit


