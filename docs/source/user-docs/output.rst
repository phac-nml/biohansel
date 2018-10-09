Output (Errors Encountered)
===========================

"WARNING - Intermediate subtype"
--------------------------------
   This warning can be triggered from four different type of conditions:
   
   **1st condition:** If more than 5% of the tiles are missing (or that less that 95% of the schemes targets are not matched) then the warning will be displayed and it will be an *error type 1: missing targets*

   **2nd condition:** If there is a clash, meaning that there are both "+" and "-" tiles present in the same genome location (above background noise level). The warning will be displayed and it will be considered an *error type 2: mixed sample*
   
   **3rd condition:** If only a fraction of the tiles are positive for the final subtype ("# of tiles matching subtype expected > # of tiles matching subtype") 
   
   **4th condition:** A warning will be produced if the tiles for the final subtype are a mixture of both "+" and "-" BUT do NOT clash for the same positions.
