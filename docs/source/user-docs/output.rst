===========================
Output (Errors Encountered)
===========================

Three different result files will be produced: tech results.tab, match results.tab & results.tab
   
  
Tech results.tab
################
================ ================================== ================================== ==================== ===========================
    `Sample`_               `Subtype`_                    `Average Tile Coverage`_        `QC status`_            `QC message`_
---------------- ---------------------------------- ---------------------------------- -------------------- ---------------------------
  (Sample Name)    (Corresponding Subtypes Found)    (Corresponding avg tile coverage) (PASS/FAIL/WARNING)  (Corresponding QC message)   
================ ================================== ================================== ==================== ===========================

Sample
------
Names of samples that are ran on biohansel

Subtype
-------
Consensus result from the analysis 
-> If mixed results, biohansel will list all different subtypes detected
-> If no "+" target is detected it will produce: "No Subtype!"

Average Tile Coverage
---------------------
This section displays the average coverage of all the targets that were present in the sample.

QC status
---------
Three possibilities based on the QC analysis described below: (hyperlink QCmsg)

PASS
FAIL
WARNING

QC message
-----------

"WARNING: Intermediate Subtype"
--------------------------------
   Warnings will be triggered if all four following conditions are met:
   
   **1st condition:** Less than 5% of the tiles are missing (or more than 95% of the schemes targets are matched)

   **2nd condition:** There should be no clash for "+" and "-" targets for the same genome position (above background noise level)
   
   **3rd condition:** Only a fraction of the tiles are positive for the final subtype ("# of tiles matching subtype expected > # of tiles matching subtype") 
   
   **4th condition:** The targets for the final subtype are a mixture of both "+" and "-" BUT do NOT clash for the same positions.

"WARNING: Low Coverage"
------------------------
If the "Avg Tile Coverage" is below the parameters given for low coverage (parameters are adjustable) (default min average coverage: 20- fold)

Error Type 1: Missing Tiles
---------------------------
*** The "Maximum amount of missing tiles to be allowed before being considered an error" can be edited based on preference and scheme
