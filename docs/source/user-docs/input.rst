Input
=====

Analysis of a single FASTA file
-------------------------------

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab /path/to/SRR1002850.fasta


Contents of ``results.tab``:

.. code-block::

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR1002850  heidelberg  2.2.2.2.1.4 2; 2.2; 2.2.2; 2.2.2.2; 2.2.2.2.1; 2.2.2.2.1.4  1037658-2.2.2.2.1.4; 2154958-2.2.2.2.1.4; 3785187-2.2.2.2.1.4   True        202 202 17  17  3   3   SRR1002850.fasta


Contents of ``match_results.tab``:

.. code-block::

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

.. code-block:: bash

    hansel -s heidelberg -vv -t 4 -o results.tab -O match_results.tab -p SRR5646583_forward.fastqsanger SRR5646583_reverse.fastqsanger


Contents of ``results.tab``:

.. code-block::

    sample  scheme  subtype all_subtypes    tiles_matching_subtype  are_subtypes_consistent inconsistent_subtypes   n_tiles_matching_all    n_tiles_matching_all_total  n_tiles_matching_positive   n_tiles_matching_positive_total n_tiles_matching_subtype    n_tiles_matching_subtype_total  file_path
    SRR5646583  heidelberg  2.2.1.1.1.1 2; 2.2; 2.2.1; 2.2.1.1; 2.2.1.1.1; 2.2.1.1.1.1  1983064-2.2.1.1.1.1; 4211912-2.2.1.1.1.1    True        202 202 20  20  2   2   SRR5646583_forward.fastqsanger; SRR5646583_reverse.fastqsanger


Contents of ``match_results.tab``:

.. code-block::

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

.. code-block:: bash

    hansel -s heidelberg -vv --threads <n_cpu> -o results.tab -O match_results.tab -D /path/to/fastas_or_fastqs/


``hansel`` will only attempt to analyze the FASTA/FASTQ files within the specified directory and will not descend into any subdirectories!
