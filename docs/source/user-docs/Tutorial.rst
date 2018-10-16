Tutorial
========

.. |heidelberg| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Specs%20for%20biohansel.PNG
   :alt: specs of biohansel run
   :width: 500 px
 
.. |experimental| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Biohansel%20location.PNG
   :alt: location of biohansel in galaxy
   :width: 250 px

**To see the different types of outputs that are produced by biohansel go to: `Output <https://bio-hansel.readthedocs.io/en/readthedocs/user-docs/output.html>`_

Testing:
########

To verify that Biohansel is working properly with the different interfaces:

Download either the CP012921.fasta (fasta file) or the SRR2598330(fastq-dump).fastqsanger.gz (raw file) which are two of the same samples:

   <https://share.corefacility.ca/index.php/s/dRGOuqhDJUNeKmE> (password: biohansel)
   
|
CP012921.fasta (fasta file) - This is the assembly file 
SRR2598330(fastq-dump).fastqsanger.gz - This is the raw reads file
\*** Both files are Salmonella heidelberg samples

**Testing Results**
-------------------

For CP012921.fasta (fasta file):
 
|
NML - Galaxy Access (BioHansel)
###############################
1.) Create a new history in Galaxy and either the fasta file or the raw reads file onto the new history
  
2.) Find Biohansel on the right-hand side in the "Tools" Section: Under the Experimental Section

  |experimental|
  
3.) For the "SNP Subtyping Scheme", select the proper scheme corresponding to the organism in your samples
       - (For verification/testing select the "Salmonella Heidelberg subtype scheme")
       
|heidelberg|
  
4.) Execute the file and three results should be produced: tech_results.tab, match_results.tab and results.tab
(If running the "testing" fasta or raw file; to verify go to `Testing results`_)

-> The .tab files can be opened in excel
