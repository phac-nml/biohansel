======
Output 
======

This page describes the three different result files will be produced from running BioHansel: `tech results.tab`_, `match results.tab`_ & `results.tab`_. The results found in these three files will be the same whether you are using the command line or Galaxy to run an analysis.


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

.. |all_subtypes| image:: all_subtypes.png
   :alt: Output of all subtypes
   :width: 450 px

.. |inconsistent_subtypes_false| image:: inconsistent_subtypes_false.png
   :alt: Output of all subtypes
   :width: 477 px

.. |matching_all| image:: matching_all.png
   :alt: tiles matching all output
   :width: 420 px

.. |good_tech| image:: good_tech.png
   :alt: Correct tech_results.tab file
   :width: 400 px

|
**Tech results.tab**
####################

Tech_results.tab is the simplest output file released by running a BioHansel analysis. It contains only the sample name, subtype, and the QC status of the sample allowing this file to be easy to interpret at the cost of not elaborating on any of the specific details of the analysis. Found below are the columns and explanations of the columns for this output file:


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
This column gives the subtype of the sample determined by the analysis. This column can display a single positive subtype, a list of positive subtypes, or no subtype depending on the results of the analysis. A good analysis will output the following:

|good_tech|

If this column does not display a single positive subtype, it will show one of the two following situations:

1. Different subtypes if mixed samples are run or there is an error in a user-created scheme. In this case, BioHansel will list all different subtypes detected.

|mixed_result|

2. If no positive target is detected, the column will be blank and the qc_message will state that no tiles/targets were found.

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

The following is the scheme for the match_results.tab file **For a single Fasta file**. **Running raw reads data has slightly different output columns due to the different nature of the data**. The output columns for the match_results.tab file are shown below broken into different charts to allow them to fit mostly on one page. In the real generated file, they would all found in the same long row. Below, you will find detailed information for each column.

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

Running raw reads files/FastQ files gives slightly different output columns when compared to the Fasta file match_results.tab output due to the slight differences in the data that each file contains. The overall output for a match_results.tab results output from a FastQ file looks as such:


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

Overall the match_results.tab file for analyzing raw reads will look as such:

+------------------------+--------------------------------+---------------------------------+-------------------------------+-------------------------+--------------+-------------------+-----------------+---------------+---------------+------------------+---------------------+----------------------------+  
| Tilename               | Sequence                       | Frequency                       | Refposition                   | Subtype                 | is_pos_tile  | is_kmer_freq_okay | File_path       | Sample        |Scheme         | Scheme_version   | QC_Status           | QC_message                 |
+========================+================================+=================================+===============================+=========================+==============+===================+=================+===============+===============+==================+=====================+============================+  
| (Name of Target/K-mer) | (Corresponding K-mer Sequence) | (Number of exact matches found) | (Match Position in reference) | (Subtypes in Tilename)  | (TRUE/FALSE) | (TRUE/FALSE)      | (File Location) | (Sample Name) |(Scheme Name)  | (Scheme Version) | (PASS/FAIL/WARNING) | (Corresponding QC message) |
+------------------------+--------------------------------+---------------------------------+-------------------------------+-------------------------+--------------+-------------------+-----------------+---------------+---------------+------------------+---------------------+----------------------------+



**Detailed Column Information** 
-------------------------------

The detailed information on the meaning of each columns outputs for both files can be found below:

Tilename
""""""""

This column gives the name of the target/kmer that matched to the sample. It will match to the name of the tile in the fasta file following the fasta convention as seen in the `input section <input.html>`_. The tiles give the identity of the sample

Sequence
""""""""

The column contains the sequence of the tile from the Tilename column. This sequence is the 33 bp fragment that matched somewhere in the sample.

is_revcomp
""""""""""

Is the tile found in the forward direction or the reverse direction?

1. FALSE - the target tile was found from the 5' to 3' direction 

2. TRUE - the target tile was found in the 3' to 5' direction in the sample

Contig_id
"""""""""

Displays the name of the contig as found in the Fasta file.

Frequency
"""""""""

Displays the exact number of matches found for the tile/k-mer in the raw reads/FastQ file input.

Match_index
"""""""""""

Displays the last nucleotide match of a k-mer/tile as its position in the genome.

For example, if the tile matched the genome from positions 12312 to 12345, the SNP would be at position 12329 and output of this column would be 12345.

Refposition
"""""""""""

Displays the numerical position of the tile/k-mers SNP in the reference genome. This information is also found in the description of the tile in the subtyping schemes Fasta file. 

Subtype
"""""""

Shows the consensus subtype of the sample as determined by the analysis. 

This column can display a single positive subtype, a list of positive subtypes, or no subtype depending on the results.

is_pos_tile
"""""""""""
Is the tile in question a positive k-mer/target for specific subtype?

1. TRUE - the positive SNP has been found in the sample

2. FALSE - the negative SNP has been found in the sample


is_kmer_freq_okay
"""""""""""""""""

Is the frequencey of the k-mer/tile within the specified QC parameters (min/max)? For FastQ datasets. 

1. TRUE - enough of the k-mer has been found in the dataset as specified by the QC parameters

2. FALSE - not enough of the k-mer has been found in the dataset as specified by the QC parameters


File path
"""""""""

The location of the input data file.


Scheme
""""""
The name of the chosen Scheme used in the analysis.

Scheme_vers
"""""""""""

The version of the chosen scheme used in the analysis.

QC Columns
""""""""""

QC Status and QC message are found in full details under their own section as they are a part of all 3 results files. This detailed information is found in the `Quality_Control`_ section.

|
**Results.tab**
################

The results.tab output file is almost exactly the same for all inputs. This file contains the overall information of the analysis and gives the final results of a BioHansel run in more detail then the tech_results.tab file. The expanded version of all information that can be obtained from this file is as such:

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

================================ ==================== ===========================
        avg_tile_coverage             QC status               QC message  
-------------------------------- -------------------- ---------------------------
(Average frequency of all tiles) (PASS/FAIL/WARNING)  (Corresponding QC message) 
================================ ==================== ===========================

Sample
------

Provides the names of samples that were run on BioHansel


Scheme
------

The name of the chosen Scheme used in the analysis.


Scheme_Version
--------------

The version of the chosen scheme used in the analysis.


Subtype
-------

Shows the consensus subtype of the sample as determined by the analysis.

This column can display a single positive subtype, a list of positive subtypes, or no subtype depending on the results.


All_subtypes
------------

All of the subtypes in all the levels of lineage leading to the final subtype.

|all_subtypes|


tiles_matching_subtype
----------------------

Displays the subtype(s) that the most downstream, specific tiles have matched to. For good, non-mixed results it should be the same as the subtype column.


are_subtypes_consistent
-----------------------

1. TRUE - the subtypes are consistent as defined.

- Consistency -> All positive tiles within QC parameters have consistent subtypes in downstream sublineages corresponding to parent subtype.

|consistent|

Each tile must become more specific to the final subtype while matching all of the previous ones to be considered consistent.

2. FALSE - the subtypes are not consistent.


inconsistent_subtypes
---------------------

If "are_subtypes_consistent" is FALSE, it lists subtypes that are inconsistent to parent.

|inconsistent_subtypes_false|


n_tiles_matching_all
--------------------

Counting all of the actual k-mer matches (both positive and negative) that make up each subtype lineage as defined by the subtyping scheme used/created.

|n_all|


n_tiles_matching_all_expected
-----------------------------

The total number k-mer/target matches expected (both positive and negative) that make up each subtype lineage as defined by the subtyping scheme used/created.

Every/almost every k-mer defined in the scheme should match somewhere in the sample if the sample is of high quality.

|matching_all|


n_tiles_matching_positive
-------------------------

The number of positive matches in the sample from all of the upstream lineages of the output subtype as defined by the subtyping scheme.

|positive|


n_tiles_matching_positive_expected
----------------------------------

The expected number of positive matches from all of the upstream lineages of the output subtype as defined by the subtyping scheme.

For a good analysis, this value should match the sample.


n_tiles_matching_subtype
------------------------

The number of positive matches in the sample sublineage only.

|subtype|


n_tiles_matching_subtype_expected
---------------------------------

The expected number of positive matches in the sample sublineage only.

File Path
---------

The file location of the input data.


Avg_tile_coverage
-----------------

The average frequency of all tiles, both positive and negative, that were found in the sample. This output column is only found for analysis of raw reads FastQ files and it is an indicator that there was a sufficient amount of overlap in the dataset for the results to be significant. 


QC Columns
----------

QC Status and QC message are found in full details under their own section as they are a part of all 3 results files. This detailed information is found in the `Quality_Control`_ section.


**Quality_Control**
###################

|
**QC Status**
-------------
Three possibilities can be shown in this column based on the QC analysis described below: `QC message`_

1. PASS

2. FAIL

3. WARNING

|
**QC message**
--------------

The QC message displayed provides information on what happened in the analysis and where, if there was a warning or fail, the data can be cleaned up/improved to obtain a passing result. 

*"Pass"*
"""""""""
A pass occurs when there is no errors in the targeted lineage and its corresponding sublineages:

|pass|

Once the QC module is declared as a pass, there is no information in the QC message column displayed. The result should be considered a valid analysis.

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
\*** The Maximum amount of missing tiles, either positive or negative, to be allowed before being considered an error/fail. This amount can be edited based on preference and scheme.

Three possible causes:

1. Bacterial scheme does not match target                                       

2. Low genome coverage or low quality data

3. Range of target coverage extends outside of QC limits (k-mer frequency thresholds default = min:8, max:500)

** To determine which cause, the average coverage depth is returned to the user. The value is calculated based on the coverage for all tiles that were above the minumum coverage threshold (indicated by the QC parameters: default value = 8) 

|missing|

|
*Error Type 2: Mixed Sample*
""""""""""""""""""""""""""""
A mixed sample error is where BioHansel is unsure what the final subtype is of the sample due to one of two possible causes:

1. BioHansel came out with an "inconsistent result" designation

2. Position conflict: both "+" and "-" targets are found in the same target genome position above background noise level

A possible solution to this error if the average genome coverage is above 100 is to increase the minimum k-mer threshold to at least 10% of the average genome coverage. This will change the background noise tolerated and potentially allow for a positive result to occur. 

|mixed|

|
*"Error Type 3: Ambiguous result"* 
""""""""""""""""""""""""""""""""""
Caused by both conditions met:

1. Total matching tiles is within 5% of the expected value

2. 3 or more tiles are missing for the final subtype call (Error 3a)

|inconsistent|

|
*"Error Type 4: Unconfident/Not confident result"*
""""""""""""""""""""""""""""""""""""""""""""""""""
Lineage call is uncertain due to missing targets in downstream sublineage.

|unconfident|

.. _schemes: subtyping_schemes.html


