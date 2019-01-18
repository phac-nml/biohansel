Input
=====

This section describes the main input files that are needed to run BioHansel along with the various methods avaliavble to analyze datasets using the tool.

The three input files are:

- The WGS data/assembled genome added as a FASTQ or FASTA file

- A chosen genotyping scheme (heidelberg or enteritidis) or a user created custom genotyping scheme (FASTA). More info in the `subtyping schemes section <subtyping_schemes.html>`_.

- An optional metadata table to add additional information to results

More detailed information on the output of what each results file contains can be found in the `Output section <output.html>`_.

Analysis of a single FASTA file
-------------------------------

Analysis of a single FASTA file would be run on a sample that has already been assembled into contigs through another program/tool. To run the command, the arguments that must be specified include:

- -s "scheme" where the scheme defined can be one of the two built or a user created one (FASTA format)

- Any combination of the results delimiters and names (file names can be changed but must be included after the argument):

    - -o results.tab
    - -O match_results.tab
    - -S tech_results.tab

- The name of the FASTA file at the end of the command

An example command for the analysis of a single FASTA file called SRR1002850.fasta would look like:

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab /path/to/SRR1002850.fasta

Or, if you have already changed to the directory containing the dataset, you can use the following command where you do not have to specifiy the path to the data:

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab SRR1002850.fasta

The output of the BioHansel tool can be found in the directory that the command was run from.

Contents of ``results.tab``:

.. code-block:: bash

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR1002850  heidelberg  2.2.2.2.1.4 2; 2.2; 2.2.2; 2.2.2.2; 2.2.2.2.1; 2.2.2.2.1.4  1037658-2.2.2.2.1.4; 2154958-2.2.2.2.1.4; 3785187-2.2.2.2.1.4   True        202 202 17  17  3   3   SRR1002850.fasta


Contents of ``match_results.tab``:

.. code-block:: bash

    tilename    stitle  pident  length  mismatch    gapopen qstart  qend    sstart  send    evalue  bitscore    qlen    slen    seq coverage    is_trunc    refposition subtype is_pos_tile sample  file_path   scheme
    775920-2.2.2.2  NODE_2_length_512016_cov_46.4737_ID_3   100.0   33  0   0   1   33  474875  474907  2.0000000000000002e-11  62.1    33  512016  GTTCAGGTGCTACCGAGGATCGTTTTTGGTGCG   1.0 False   775920  2.2.2.2 True    SRR1002850  SRR1002850.fasta   heidelberg
    negative3305400-2.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  276235  276267  2.0000000000000002e-11  62.1    33  427905  CATCGTGAAGCAGAACAGACGCGCATTCTTGCT   1.0 False   negative3305400 2.1.1.1 False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3200083-2.1 NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  170918  170950  2.0000000000000002e-11  62.1    33  427905  ACCCGGTCTACCGCAAAATGGAAAGCGATATGC   1.0 False   negative3200083 2.1 False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3204925-2.2.3.1.5   NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  175760  175792  2.0000000000000002e-11  62.1    33  427905  CTCGCTGGCAAGCAGTGCGGGTACTATCGGCGG   1.0 False   negative3204925 2.2.3.1.5   False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3230678-2.2.2.1.1.1 NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  201513  201545  2.0000000000000002e-11  62.1    33  427905  AGCGGTGCGCCAAACCACCCGGAATGATGAGTG   1.0 False   negative3230678 2.2.2.1.1.1 False   SRR1002850  SRR1002850.fasta   heidelberg
    negative3233869-2.1.1.1.1   NODE_3_length_427905_cov_48.1477_ID_5   100.0   33  0   0   1   33  204704  204736  2.0000000000000002e-11  62.1    33  427905  CAGCGCTGGTATGTGGCTGCACCATCGTCATTA   1.0 False   
    [Next 196 lines omitted.]


Analysis of a single FASTQ readset
----------------------------------

Analysis of a single FASTQ readset would be run on raw sequencing data. To run the command, the arguments that must be specified include:

- -s "scheme" where the scheme defined can be one of the two built or a user created one (FASTA format)

- Any combination of the results delimiters and names (file names can be changed but must be included after the argument):
 
    - -o results.tab
    - -O match_results.tab
    - -S tech_results.tab

- The name of the FASTQ file(s) at the end of the command
    - For single-end reads include the one file 
    - For paired-end reads include both files one after the other

An example command for the analysis of a single paired-end reads run dataset would look like:

.. code-block:: bash

    hansel -s heidelberg -vv -t 4 -o results.tab -O match_results.tab -p SRR5646583_forward.fastqsanger SRR5646583_reverse.fastqsanger


Contents of ``results.tab``:

.. code-block:: bash

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR5646583  heidelberg  2.2.1.1.1.1 2; 2.2; 2.2.1; 2.2.1.1; 2.2.1.1.1; 2.2.1.1.1.1  1983064-2.2.1.1.1.1; 4211912-2.2.1.1.1.1    True     single   202 202 20  20  2   2   SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger


Contents of ``match_results.tab``:

.. code-block:: bash

    seq freq    sample  file_path   tilename    is_pos_tile subtype refposition is_kmer_freq_okay   scheme
    ACGGTAAAAGAGGACTTGACTGGCGCGATTTGC   68  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    21097-2.2.1.1.1 True    2.2.1.1.1   21097   True    heidelberg
    AACCGGCGGTATTGGCTGCGGTAAAAGTACCGT   77  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    157792-2.2.1.1.1    True    2.2.1.1.1   157792  True    heidelberg
    CCGCTGCTTTCTGAAATCGCGCGTCGTTTCAAC   67  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    293728-2.2.1.1  True    2.2.1.1 293728  True    heidelberg
    GAATAACAGCAAAGTGATCATGATGCCGCTGGA   91  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    607438-2.2.1    True    2.2.1   607438  True    heidelberg
    CAGTTTTACATCCTGCGAAATGCGCAGCGTCAA   87  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    691203-2.2.1.1  True    2.2.1.1 691203  True    heidelberg
    CAGGAGAAAGGATGCCAGGGTCAACACGTAAAC   33  SRR5646583 SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger    944885-2.2.1.1.1    True    2.2.1.1.1   944885  True    heidelberg
    [Next 200 lines omitted.]

Analysis of all FASTA/FASTQ files in a directory
------------------------------------------------

Analysis on **all** of the FASTA/FASTQ files in the specified directory. This will run on all FASTA/FASTQ files in the directory so be sure that there are no miscellaneous files that may increase analysis time or lead to errors.

``bio_hansel`` will only attempt to analyze the FASTA/FASTQ files within the specified directory and will not descend into any subdirectories! As such, make sure all of the data to be analyzed is in the same location.

Analysis of all of the sequencing files in a directory must include following the arguments to run properly:

- -s "scheme" where the scheme defined can be one of the two built or a user created one (FASTA format)

- --threads <#_cpu> to specify the number of CPU's wanted to run the analysis. (give a number)

- Any combination of the results delimiters and names (file names can be changed but must be included after the argument):
 
    - -o results.tab
    - -O match_results.tab
    - -S tech_results.tab

- -D /path/to/directory_with_data

An example of a general command for the analysis of a directory of FASTA/FASTQ files:

.. code-block:: bash

    hansel -s heidelberg -vv --threads <n_cpu> -o results.tab -O match_results.tab -D /path/to/fastas_or_fastqs/

The chosen output files can still be found in the directory that the command was run from and will contain data from each of the analyzed files run by bio_hanzel. 

Ex. If you had your data directory in the path USER/name of user/bio_hansel/data and ran your command in the USER/name of user/bio_hansel folder, then the results of the analysis would end up in the bio_hansel folder. To run the analysis on the folder in this situation, your command would be as follows:

.. code-block:: bash

    hansel -s heidelberg -vv --threads 1 -o results.tab -O match_results.tab -D data/


Subtype Metadata Table (Optional)
---------------------------------

Optionally you can select a subtype metadata information table to include subtype metadata along with the subtyping results created with BioHansel. Metadata tables must be in a tab-delimited format to correctly work. The file extension for your metadata table should be **.tsv** if at all possible or you may end up with an error.

To add a metadata table to the analysis you will add the argument `-M <metadata_scheme.tsv>` to any other analysis command. There are no requirements for the number of columns or the content of each of the columns on the metadata table so long as the first column is labeled as "subtype". 

A command that incorporates the -M command for analysis would be structured following the previously established requirements and looks as follows:

.. code-block:: bash

    hansel -s heidelberg -M <metadata_scheme.tsv> -vv -o results.tab -O match_results.tab <data>


The BioHansel results table will be joined with the subtype metadata table based if a subtype on the metadata table matches one on the results. If a match occurs, the metadata of that subtype will be added to the table at the end of the results.tab and tech_results.tab results files. 

Example metadata table (called meta.tsv):

+-------------+-------+--------+----------+ 
| Subtype     | Clade | Source | Symptoms | 
+=============+=======+========+==========+  
| 1           | I     | Geese  | Death    | 
+-------------+-------+--------+----------+ 
| 1.1         | I     | Moose  | Burns    | 
+-------------+-------+--------+----------+ 
| 2.2.1.1.1   | II    | Mouse  | Boils    | 
+-------------+-------+--------+----------+  
| 2.2.2.2.2.1 | IIa   | Human  | Rash     | 
+-------------+-------+--------+----------+ 

***When naming a metadata table make sure there are no spaces or parentheses and that its extension is .tsv or the analysis may fail.*** 

The added metadata will appear at then end of the results.tab and the tech_results.tab files.

Example: tech_results.tab without metadata added:

+----------+-----------+-----------+------------+ 
| Sample   | Subtype   | qc_status | qc_message | 
+==========+===========+===========+============+  
| CP012921 | 2.2.3.1.2 | PASS      |            | 
+----------+-----------+-----------+------------+ 

tech_results.tab with metadata:

+----------+-----------+-----------+------------+-------+--------+----------+ 
| Sample   | Subtype   | qc_status | qc_message | Clade | Source | Symptoms | 
+==========+===========+===========+============+=======+========+==========+  
| CP012921 | 2.2.3.1.2 | PASS      |            | I     | Geese  | Rash     | 
+----------+-----------+-----------+------------+-------+--------+----------+ 

You can add metadata to the analysis with Galaxy by uploading either a .tsv or a .csv file to your history and specifying that you want it used in the analysis. A .tsv file is recommended.






