======
Output 
======

This page describes the three different result files will be produced from running BioHansel: `tech results.tab`_, `match results.tab`_ & `results.tab`_. These results will be the same whether you are using the command line or Galaxy to run an analysis.


.. |mixed| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Mixed.PNG
   :width: 100 px
   :alt: Example of Mixed targets
   
   
.. |missing| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Missing%20Targets.PNG
   :width: 100 px
   :alt: Example of Missing Targets
   
.. |inconsistent| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Inconsistent%20results.PNG
   :width: 100 px
   :alt: Example of Inconsistent results
   
.. |unconfident| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Unconfident%20(1).PNG
   :width: 100 px
   :alt: Example of Unconfident results
   
.. |pass| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Pass.PNG
   :alt: This is an ideal picture of a passed scheme
   :width: 100 px
Incorrect

.. |positive| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Positive%20pic%20of%20matching.PNG
   :alt: picture of positive match
   :width: 100 px

.. |consistent| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/PCIS%20BIO.PNG
   :alt: picture of consistent
   :width: 100 px

.. |n_all| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/N%20tiles%20all%20picture.PNG
   :alt: picture of all match
   :width: 100 px
 
.. |subtype| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/sUBTYPE%20MATCH%20PIC.PNG
   :alt: picture of subtype match
   :width: 100 px

.. |mixed_result| image:: mixed_sub_result.png
   :alt: Mixed subtype result
   :width: 500 px

.. |error_no_result| image:: No_result.png
   :alt: no result
   :width: 600 px


|
**Tech results.tab**
####################

Tech results is the simplest output file and only contains the sample name, the subtype found by the analysis, and the quality control status with message. Here the output file :


+---------------+--------------------------------+-------------------------------------------+
| Sample        | Subtype                        | Avg_tile_coverage                         |
+===============+================================+===========================================+
| (Sample Name) | (Corresponding Subtypes Found) | (average tile coverage of all the targets |
+---------------+--------------------------------+-------------------------------------------+

+---------------------+----------------------------+ 
| QC_status           | QC_message                 |
+=====================+============================+ 
| (PASS/FAIL/WARNING) | (Corresponding QC message) |
+---------------------+----------------------------+


Sample
------
This column provides the names of samples that were run on BioHansel


Subtype
-------
This column gives the subtype of the sample determined by the analysis. This column can display a single positive subtype, a list of positive subtypes, or no subtype depending on the results of the analysis. If this column does not display a single positive subtype, it will show one of the two following situations:

- Different subtypes if mixed samples are run or there is an error in a user-created scheme. In this case, BioHansel will list all different subtypes detected.

|mixed_result|

- If no positive target is detected, the column will be blank and the qc_message will state that no tiles/targets were found.

|error_no_result|

Average Tile Coverage
---------------------

Found only when analyzing raw read FastQ files. It displays the average coverage of all of the targets/k-mers that were present in the sample.

QC Columns
----------

QC Status and QC message are found in full details under their own section as they are a part of all 3 results files. This detailed information is found in the `Quality_Control`_ section.



**Match Results.tab**
#####################


**Fasta File Output**
---------------------

The following is the scheme for the match_results.tab file **For a single Fasta file**. **Running raw reads data has slightly different output columns**. The columns for a contig contained in a are broken up to fit all on one page with out scrolling. Here, you can see all of the outputs of this file at a glance along with a short block of info about the column.

+------------------------+--------------------------------+--------------+------------------+ 
| Tilename               | Sequence                       | is_revcomp   | Contig_id        |
+========================+================================+==============+==================+  
| (Name of Target/K-mer) | (Corresponding K-mer Sequence) | (TRUE/FALSE) | (Name of Contig) |
+------------------------+--------------------------------+--------------+------------------+


+------------------+-------------------------------+-------------------------+--------------+ 
| Match_index      | Refposition                   | Subtype                 | is_pos_tile  |
+==================+===============================+=========================+==============+  
| (Match Position) | (Match Position in reference) | (Subtypes in Tilename)  | (TRUE/FALSE) |
+------------------+-------------------------------+-------------------------+--------------+


+---------------+-----------------+---------------+------------------+
| Sample        | File_path       | Scheme        | Scheme_version   |
+===============+=================+===============+==================+ 
| (Sample Name) | (File Location) | (Scheme Name) | (Scheme Version) |
+---------------+-----------------+---------------+------------------+


+---------------------+----------------------------+ 
| QC_status           | QC_message                 |
+=====================+============================+ 
| (PASS/FAIL/WARNING) | (Corresponding QC message) |
+---------------------+----------------------------+

All of the columns in the correct order in the match_results.tab file looks as such:

+------------------------+--------------------------------+--------------+------------------+------------------+-------------------------------+-------------------------+--------------+---------------+-----------------+---------------+------------------+---------------------+----------------------------+  
| Tilename               | Sequence                       | is_revcomp   | Contig_id        | Match_index      | Refposition                   | Subtype                 | is_pos_tile  | Sample        | File_path       | Scheme        | Scheme_version   | QC_Status           | QC_message                 |
+========================+================================+==============+==================+==================+===============================+=========================+==============+===============+=================+===============+==================+=====================+============================+ 
| (Name of Target/K-mer) | (Corresponding K-mer Sequence) | (TRUE/FALSE) | (Name of Contig) | (Match Position) | (Match Position in reference) | (Subtypes in Tilename)  | (TRUE/FALSE) | (Sample Name) | (File Location) | (Scheme Name) | (Scheme Version) | (PASS/FAIL/WARNING) | (Corresponding QC message) |
+------------------------+--------------------------------+--------------+------------------+------------------+-------------------------------+-------------------------+--------------+---------------+-----------------+---------------+------------------+---------------------+----------------------------+

**Raw Reads FastQ File Output**
-------------------------------

The raw reads file gives slightly different outputs when compared to the Fasta file match_results.tab output. The overall Output looks as such:


+------------------------+--------------------------------+---------------------------------+ 
| Tilename               | Sequence                       | Frequency                       |
+========================+================================+=================================+  
| (Name of Target/K-mer) | (Corresponding K-mer Sequence) | (Number of exact matches found) |
+------------------------+--------------------------------+---------------------------------+


+-------------------------------+-------------------------+--------------+-------------------+ 
| Refposition                   | Subtype                 | is_pos_tile  | is_kmer_freq_okay |
+===============================+=========================+==============+===================+  
| (Match Position in reference) | (Subtypes in Tilename)  | (TRUE/FALSE) | (TRUE/FALSE)      |
+-------------------------------+-------------------------+--------------+-------------------+


+-----------------+---------------+---------------+------------------+
| File_path       | Sample        |Scheme         | Scheme_version   |
+=================+===============+===============+==================+ 
| (File Location) | (Sample Name) |(Scheme Name)  | (Scheme Version) |
+-----------------+---------------+---------------+------------------+


+---------------------+----------------------------+ 
| QC_Status           | QC_message                 |
+=====================+============================+ 
| (PASS/FAIL/WARNING) | (Corresponding QC message) |
+---------------------+----------------------------+

**Detailed Column Information** 
-------------------------------

Tilename
""""""""

This column gives the name of the target/kmer that matched to the sample. It will match to the name of the tile in the fasta file following the fasta convention as seen in the `input section <input.html>`_. The tiles give the identity of the sample

Sequence
""""""""

This column gives the sequence of the tile from the Tilename column. This sequence is the 33 bp fragment that matched somewhere in the sample.

is_revcomp
""""""""""

If the target tile was found from the 5' to 3' direction in the sample, this column will display "FALSE".

If the target tile was found in the 3' to 5' direction in the sample, this column will display "TRUE".

Contig_id
"""""""""

This column displays the name of the contig as found in the Fasta file.

Frequency
"""""""""



Match_index
"""""""""""

This column outputs the last nucleotide match of a k-mer/tile as its position in the genome.

For example, if the tile matched the genome from positions 12312 to 12345, the SNP would be at position 12329 and output of this column would be 12345.

Refposition
"""""""""""

This column displays the reference position of the SNP as found in the k-mer tile.

Subtype
"""""""

This column gives the subtype of the sample determined by the analysis. This column can display a single positive subtype, a list of positive subtypes, or no subtype depending on the results of the analysis. If this column does not display a single positive subtype, it will show one of the two following situations:

- Different subtypes if mixed samples are run or there is an error in a user-created scheme. In this case, BioHansel will list all different subtypes detected.


is_pos_tile
"""""""""""
Is it a positive k-mer/target for specific subtype?

1.) TRUE

2.) FALSE


is_kmer_freq_okay
"""""""""""""""""

Is it within the specified QC parameters (min/max)

1.) TRUE

2.) FALSE


File path
"""""""""

File Location


Scheme
""""""
Name of the given Scheme

Scheme_vers
"""""""""""

Version of the given scheme

|
**Results.tab**
################

===================== ======================= =============================== ========================== ============================
       Sample                Sequence                  Scheme_vers                    Subtype                  all_subtype  
--------------------- ----------------------- ------------------------------- -------------------------- ----------------------------
    (Sample Name)          (Scheme name)            (Version of Scheme)         (Subtypes in tilename)    (Subtypes in all lineages)
===================== ======================= =============================== ========================== ============================

==================================== ============================== =========================== =======================================
     tiles_matching_subtype             are_subtypes_consistent        inconsistent_subtypes              n_tiles_matching_all   
------------------------------------ ------------------------------ --------------------------- ---------------------------------------
 (subtypes that match given tiles)            (TRUE/FALSE)                  (TRUE/FALSE)          (Number of actual matches in sample)
==================================== ============================== =========================== =======================================
 
====================================== ========================================= ========================================
    n_tiles_matching_all_expected               n_tiles_matching_positive           n_tiles_matching_positive_expected       
-------------------------------------- ----------------------------------------- ----------------------------------------
(Expected positive matches in sample)   (Number of matches in targeted lineage)   (Expected matches in targeted lineage)          
====================================== ========================================= ========================================

============================================ =========================================== =====================
      n_tiles_matching_subtype                    n_tiles_matching_subtype_expected           File path   
-------------------------------------------- ------------------------------------------- ---------------------
(Number of matches in specific sublineage)    (Expected matches in targeted sublineage)    (File Location)         
============================================ =========================================== =====================

==================== ===========================
     QC status               QC message  
-------------------- ---------------------------
 (PASS/FAIL/WARNING)  (Corresponding QC message) 
==================== ===========================


all_subtype
-----------
All of the subtypes in all the levels of lineage


tiles_matching_subtype
----------------------
(blank)


are_subtypes_consistent
-----------------------
- Consistency -> All positive tiles within QC parameters, have consistent subtypes in downstream sublineages corresponding to parent subtype

|consistent|

inconsistent_subtypes
---------------------
If "are_subtypes_consistent" is FALSE, it lists subtypes that are inconsistent to parent


n_tiles_matching_all
--------------------
Counting actual positive matches per subtype found in sample based on subtype scheme in all lineages

|n_all|


n_tiles_matching_all_expected
-----------------------------
The number positive matches expected per subtype found in sample based on subtype scheme


n_tiles_matching_positive
-------------------------
The number of positive matches in the full sample lineage 

|positive|


n_tiles_matching_positive_expected
----------------------------------
The number of positive matches expected in the full sample lineage 

n_tiles_matching_subtype
------------------------
The number of positive matches in the sample sublineage only

|subtype|

n_tiles_matching_subtype_expected
---------------------------------
The number of positive matches expected in the sample sublineage only


**Quality_Control**
###################

|
**QC Status**
-------------
Three possibilities based on the QC analysis described below: `QC message`_

1.) PASS

2.) FAIL

3.) WARNING

|
**QC message**
--------------

*"Pass"*
"""""""""
A pass occurs when there is no errors in the targeted lineage and its corresponding sublineages:

|pass|

|
*"WARNING: Intermediate Subtype"*
"""""""""""""""""""""""""""""""
Warnings will be triggered if all four following conditions are met:
   
**1st condition:** Less than 5% of the tiles are missing (by default) or more than 95% of the schemes targets are matched (parameters for this is adjustable prior to running biohansel)

**2nd condition:** There should be no clash for "+" and "-" targets for the same genome position (above background noise level)
   
**3rd condition:** Only a fraction of the tiles are positive for the final subtype ("# of tiles matching subtype expected > # of tiles matching subtype") 
   
**4th condition:** The targets for the final subtype are a mixture of both "+" and "-" BUT do NOT clash for the same positions.

|
*"WARNING: Low Coverage"*
"""""""""""""""""""""""
If the "Avg Tile Coverage" is below the parameters given for low coverage (parameters are adjustable) (default min average coverage: 20- fold)

Average coverage calculated from all targets found in the sample (The value is returned to the user)

|
*Error Type 1: Missing Tiles*
"""""""""""""""""""""""""""""
\*** The "Maximum amount of missing tiles to be allowed before being considered an error" can be edited based on preference and scheme

Two possible causes:

1.) Bacterial scheme does not match target                                       

2.) Low genome coverage or low quality data

3.) Range of target coverage extends outside of QC limits (k-mer frequency thresholds default = min:8, max:500)

** To determine which cause, the average coverage depth is returned to the user. The value is calculated based on the coverage for all tiles that were above the minumum coverage threshold (indicated by the QC parameters: default value = 8) 

|missing|

|
*Error Type 2: Mixed Sample*
""""""""""""""""""""""""""""
Two possible causes:

1.) BioHansel came out with an "inconsistent result" designation

2.) Position conflict: both "+" and "-" targets are found in the same target genome position above background noise level
-> (possible solution) if the average genome coverage is above 100, increase the minimum k-mer threshold to at least 10% of the average genome coverage

|mixed|

|
*"Error Type 3: Ambiguous result"* 
""""""""""""""""""""""""""""""""""
Caused by both conditions met:

1.) Total matching tiles is within 5% of the expected value
2.) 3 or more tiles are missing for the final subtype call (Error 3a)

|inconsistent|

|
*"Error Type 4: Unconfident/Not confident result"*
""""""""""""""""""""""""""""""""""""""""""""""""""
Lineage call is uncertain due to missing targets in downstream sublineage

|unconfident|

.. _schemes: subtyping_schemes.thml


