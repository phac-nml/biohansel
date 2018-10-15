======
Output 
======

Three different result files will be produced: `tech results.tab`_, `match results.tab`_ & `results.tab`_

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
|
**Tech results.tab**
###################
================ ================================== ================================== ==================== ===========================
    `Sample`_               `Subtype`_                    `Average Tile Coverage`_        `QC status`_            `QC message`_
---------------- ---------------------------------- ---------------------------------- -------------------- ---------------------------
  (Sample Name)    (Corresponding Subtypes Found)    (Corresponding avg tile coverage) (PASS/FAIL/WARNING)  (Corresponding QC message)   
================ ================================== ================================== ==================== ===========================

**Sample**
----------
Names of samples that are ran on biohansel


**Subtype**
-----------
Consensus result from the analysis 
If mixed results, biohansel will list all different subtypes detected
If no "+" target is detected it will produce: "No Subtype!"


**Average Tile Coverage**
-------------------------
This section displays the average coverage of all the targets that were present in the sample.

|
**Match Results.tab**
#####################

===================== ============================ =============================== ========================= =========================
    `Tilename`_               `Sequence`_                 `Frequency`_                  `ref_position`_            `Subtype`_
--------------------- ---------------------------- ------------------------------- ------------------------- -------------------------
(Name of target/kmer)   (Corresponding Sequence)   (Number of exact matches found) (pos in reference genome)  (Subtypes in tilename)   
===================== ============================ =============================== ========================= =========================

================== ============================== =========================== ======================== ========================
  `is_pos_tile`_         `is_kmer_freq_okay`_            `File path`_                 `Sample`_                `Scheme`_
------------------ ------------------------------ --------------------------- ------------------------ ------------------------
   (TRUE/FALSE)             (TRUE/FALSE)                (File location)             (Sample Name)            (Scheme Name)   
================== ============================== =========================== ======================== ========================

================= ===================== ============================
  `Scheme_vers`_       `QC status`_            `QC message`_
----------------- --------------------- ----------------------------
(Scheme version)   (PASS/FAIL/WARNING)   (Corresponding QC message)
================= ===================== ============================


**Tilename**
------------
Name of the target/kmer


**Sequence**
-----------
Corresponding sequence to target


**Frequency**
-------------
Number of exact matches found 


**ref_position**
---------------
Specific position in the reference genome used


**is_pos_tile**
---------------
Is it a positive k-mer/target for specific subtype?

1.) TRUE

2.) FALSE


**is_kmer_freq_okay**
---------------------
Is it within the specified QC parameters (min/max)

1.) TRUE

2.) FALSE


**File path**
-------------
File Location


**Scheme**
----------
Name of the given Scheme

**Scheme_vers**
---------------
Version of the given scheme

|
**Results.tab**
################

===================== ======================= =============================== ========================== ============================
      `Sample`_             `Sequence`_               `Scheme_vers`_                 `Subtype`_               `all_subtype`_
--------------------- ----------------------- ------------------------------- -------------------------- ----------------------------
    (Sample Name)          (Scheme name)            (Version of Scheme)         (Subtypes in tilename)    (Subtypes in all lineages)
===================== ======================= =============================== ========================== ============================

==================================== ============================== =========================== =======================================
    `tiles_matching_subtype`_          `are_subtypes_consistent`_     `inconsistent_subtypes`_           `n_tiles_matching_all`_
------------------------------------ ------------------------------ --------------------------- ---------------------------------------
 (subtypes that match given tiles)            (TRUE/FALSE)                  (TRUE/FALSE)          (Number of actual matches in sample)
==================================== ============================== =========================== =======================================
 
====================================== ========================================= ========================================
   `n_tiles_matching_all_expected`_            `n_tiles_matching_positive`_        `n_tiles_matching_positive_expected`_     
-------------------------------------- ----------------------------------------- ----------------------------------------
(Expected positive matches in sample)   (Number of matches in targeted lineage)   (Expected matches in targeted lineage)          
====================================== ========================================= ========================================

============================================ =========================================== =====================
     `n_tiles_matching_subtype`_                 `n_tiles_matching_subtype_expected`_        `File path`_ 
-------------------------------------------- ------------------------------------------- ---------------------
(Number of matches in specific sublineage)    (Expected matches in targeted sublineage)    (File Location)         
============================================ =========================================== =====================

==================== ===========================
    `QC status`_            `QC message`_
-------------------- ---------------------------
 (PASS/FAIL/WARNING)  (Corresponding QC message) 
==================== ===========================


**all_subtype**
---------------
All of the subtypes in all the levels of lineage


**tiles_matching_subtype**
--------------------------
(blank)


**are_subtypes_consistent**
---------------------------
- Consistency -> All positive tiles within QC parameters, have consistent subtypes in downstream sublineages corresponding to parent subtype

|consistent|

**inconsistent_subtypes**
-------------------------
If "are_subtypes_consistent" is FALSE, it lists subtypes that are inconsistent to parent


**n_tiles_matching_all**
-------------------------
Counting actual positive matches per subtype found in sample based on subtype scheme in all lineages

|n_all|


**n_tiles_matching_all_expected**
---------------------------------
The number positive matches expected per subtype found in sample based on subtype scheme


**n_tiles_matching_positive**
------------------------------
The number of positive matches in the full sample lineage 

|positive|


**n_tiles_matching_positive_expected**
--------------------------------------
The number of positive matches expected in the full sample lineage 

**n_tiles_matching_subtype**
-----------------------------
The number of positive matches in the sample sublineage only

|subtype|

**n_tiles_matching_subtype_expected**
-------------------------------------
The number of positive matches expected in the sample sublineage only

|
**QC status**
-------------
Three possibilities based on the QC analysis described below: `QC message`_

1.) PASS

2.) FAIL

3.) WARNING

|
**QC message**
---------------

|
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
"""""""""""""""""""""""""""
*** The "Maximum amount of missing tiles to be allowed before being considered an error" can be edited based on preference and scheme

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
