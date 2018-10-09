Output (Errors Encountered)
===========================

Table of Contents
#################
1.) `"WARNING - Intermediate subtype"`_

   2.) `Error Type 1: Missing Tiles`_
   
  
  
========== ========= ======================= ============ ============
  sample    subtype   average_tile_coverage   qc_status    qc_message
---------- --------- ----------------------- ------------ ------------
EC20120483 2.2.3.1.2          81.628              PASS         hi    
========== ========= ======================= ============ ============



Sample
------


Subtype
-------


Average Tile Coverage
---------------------

QC Status
---------

QC  Message
-----------


"WARNING - Intermediate subtype"
--------------------------------
   This warning can be triggered from four different type of conditions:
   
   **1st condition:** If more than 5% of the tiles are missing (or that less that 95% of the schemes targets are not matched) then the warning will be displayed and it will be an *error type 1: missing targets*

   **2nd condition:** If there is a clash, meaning that there are both "+" and "-" tiles present in the same genome location (above background noise level). The warning will be displayed and it will be considered an *error type 2: mixed sample*
   
   **3rd condition:** If only a fraction of the tiles are positive for the final subtype ("# of tiles matching subtype expected > # of tiles matching subtype") 
   
   **4th condition:** A warning will be produced if the tiles for the final subtype are a mixture of both "+" and "-" BUT do NOT clash for the same positions.


Error Type 1: Missing Tiles
---------------------------
*** The "Maximum amount of missing tiles to be allowed before being considered an error" can be edited based on preference and scheme
