======
Output 
======

Three different result files will be produced: tech results.tab, match results.tab & results.tab

.. |mixed| image:: https://github.com/phac-nml/biohansel/blob/readthedocs/docs/source/user-docs/Mixed.PNG
   :width: 100 px
.. |missing| image:: https://github.com/phac-nml/biohansel/blob/readthedocs/docs/source/user-docs/Missing%20Targets.PNG
   :width: 100 px
.. |inconsistent| image:: https://github.com/phac-nml/biohansel/blob/readthedocs/docs/source/user-docs/Inconsistent%20results.PNG
   :width: 100 px
.. |unconfident| image:: https://github.com/phac-nml/biohansel/blob/readthedocs/docs/source/user-docs/Unconfident.PNG
   :width: 100 px
.. |pass| image:: https://github.com/phac-nml/biohansel/blob/readthedocs/docs/source/user-docs/Pass.PNG


Tech results.tab
################
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




**QC status**
-------------
Three possibilities based on the QC analysis described below: `QC message`_

1.) PASS

2.) FAIL

3.) WARNING




**QC message**
---------------

|pass|

This is an ideal picture of a passed scheme

*"WARNING: Intermediate Subtype"*
"""""""""""""""""""""""""""""""
Warnings will be triggered if all four following conditions are met:
   
**1st condition:** Less than 5% of the tiles are missing (by default) or more than 95% of the schemes targets are matched (parameters for this is adjustable prior to running biohansel)

**2nd condition:** There should be no clash for "+" and "-" targets for the same genome position (above background noise level)
   
**3rd condition:** Only a fraction of the tiles are positive for the final subtype ("# of tiles matching subtype expected > # of tiles matching subtype") 
   
**4th condition:** The targets for the final subtype are a mixture of both "+" and "-" BUT do NOT clash for the same positions.


*"WARNING: Low Coverage"*
"""""""""""""""""""""""
If the "Avg Tile Coverage" is below the parameters given for low coverage (parameters are adjustable) (default min average coverage: 20- fold)

Average coverage calculated from all targets found in the sample (The value is returned to the user)


*Error Type 1: Missing Tiles*
"""""""""""""""""""""""""""
*** The "Maximum amount of missing tiles to be allowed before being considered an error" can be edited based on preference and scheme

Two possible causes:

1.) Bacterial scheme does not match target                                       

2.) Low genome coverage or low quality data

3.) Range of target coverage extends outside of QC limits (k-mer frequency thresholds default = min:8, max:500)

** To determine which cause, the average coverage depth is returned to the user. The value is calculated based on the coverage for all tiles that were above the minumum coverage threshold (indicated by the QC parameters: default value = 8) 

|missing|                                                                                                                                                                                                                                                                                                  

*Error Type 2: Mixed Sample*
""""""""""""""""""""""""""""
Two possible causes:

1.) BioHansel came out with an "inconsistent result" designation

2.) Position conflict: both "+" and "-" targets are found in the same target genome position above background noise level
-> (possible solution) if the average genome coverage is above 100, increase the minimum k-mer threshold to at least 10% of the average genome coverage

|mixed|



*"Error Type 3: Ambiguous result"* 
""""""""""""""""""""""""""""""""""
Caused by both conditions met:

1.) Total matching tiles is within 5% of the expected value
2.) 3 or more tiles are missing for the final subtype call (Error 3a)

|inconsistent|


*"Error Type 4: Unconfident/Not confident result"*
""""""""""""""""""""""""""""""""""""""""""""""""""
Lineage call is uncertain due to missing targets in downstream sublineage

|unconfident|
