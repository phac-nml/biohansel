Tutorial
========

.. |heidelberg| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Specs%20for%20biohansel.PNG
   :alt: specs of biohansel run
   :width: 500 px
 
.. |experimental| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Biohansel%20location.PNG
   :alt: location of biohansel in galaxy
   :width: 250 px
   
   
.. |fmatch| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Match_results.PNG
   :alt: fasta match results
   :width: 670 px
   
.. |ftech| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/tech_results.PNG
   :alt: fasta tech results
   :width: 600 px
   
.. |fresults| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Results.PNG
   :alt: fasta results
   :width: 900 px
   
   
.. |rmatch| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Match%20results.PNG
   :alt: raw match
   :width: 600 px
   
   
.. |rresults| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/results.PNG
   :alt: raw results
   :width: 600 px
   
   
.. |rtech| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Tech%20resultss.PNG
   :alt:  raw tech results
   :width: 600 px

.. |command| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Screen%20Shot%202018-10-18%20at%203.22.52%20PM.png
   :alt: command line commands
   :width: 600 px
   


\**To see the different types of outputs that are produced by biohansel go to: `Output <https://bio-hansel.readthedocs.io/en/readthedocs/user-docs/output.html>`_

Testing:
########

To verify that Biohansel is working properly with the different interfaces:

Download either the CP012921.fasta (fasta file) or the SRR2598330(fastq-dump).fastqsanger.gz (raw file) which are two of the same samples:

   <https://share.corefacility.ca/index.php/s/dRGOuqhDJUNeKmE> (password: biohansel)
   
CP012921.fasta (fasta file) - This is the assembly file 
SRR2598330(fastq-dump).fastqsanger.gz - This is the raw reads file
\*** Both files are Salmonella heidelberg samples

**Testing Results**
-------------------

**For CP012921.fasta (fasta file):**

*Fasta match_result.tab:*

|fmatch|

*Fasta tech_result.tab:*

|ftech|

*Fasta result.tab:*

|fresults|

|
**For SRR2598330(fastq-dump).fastqsanger.gz (raw file):**

*Raw/FASTQ match_result.tab:*

|rmatch|

*Raw/FASTQ tech_result.tab:*

|rtech|

*Raw/FASTQ result.tab:*

|rresults|

|
NML - Galaxy Access (BioHansel)
###############################
1.) Create a new history in Galaxy and either the fasta file or the raw reads file onto the new history
  
2.) Find Biohansel on the right-hand side in the "Tools" Section: Under the Experimental Section

  |experimental|
  
3.) For the "Sequence Data Type", select the proper type of data (FASTA vs. FASTAQ (raw))

4.) For the "SNP Subtyping Scheme", select the proper scheme corresponding to the organism in your samples

    (For verification/testing select the "Salmonella Heidelberg subtype scheme")
       
|heidelberg|
  
5.) Execute the file and three results should be produced: tech_results.tab, match_results.tab and results.tab
(If running the "testing" fasta or raw file; to verify go to `Testing results`_)

-> The .tab files can be opened in excel


Running BioHansel on Terminal (MAC) using Conda
###############################################
1.) Go to `Installation <https://bio-hansel.readthedocs.io/en/readthedocs/user-docs/usage.html>`_ and download Miniconda from the website following the instructions corresponding to your given iOS

2.) After installing Conda, go on terminal and create a conda environment by inputing this command:

conda create -n *name of environment* python=3.6

3.) It will ask you to proceed (y/n) afterwards, type in: y

4.) Then activate your environment by typing:

source activate *name of your environment*

|
5.) Now install biohansel onto conda environment by inputting:

conda install -c bioconda bio_hansel

|
6.) To confirm that biohansel has been installed in the environment, input:

biohansel -h 

#this command shows the numerous types of commands you can use in for biohansel

go to `command-line <https://bio-hansel.readthedocs.io/en/readthedocs/user-docs/command-line.html>`_ to see detailed description

|command|

|
7.) Then input:

pwd 

#pwd command stands for print working directory, which shows what directory you are currently in

|
8.) Using the directory you are in (which is most likely User/"*name of user*) you point the terminal to go to the directory where the file is by inputting:

cd *where the file is*

Example: (if the file was in User/name of user/Downloads) you input:

cd User/name of user/Downloads

# cd (change directory) command

|
9.) This will put you straight into the directory where the file is. Then just run the file using this as an example:

hansel -s heidelberg -vv -o results.tab -O match_results.tab -S tech_results.tab CP012921.fasta

-s -> this command is for the name of the scheme used in biohansel (enteritidis and heidelberg are the two built in schemes right now)

-o -> this command is for the most basic of results (you can change the name to whatever you want *just remember to add .tab)

-O -> this command is for a more detailed type of results (known as match_results.tab, but you can change it to whatever name you want)

-S -> this command is for the tech_results.tab (change name to whatever you want *just remember to add .tab)

Then at the end of the command just input the name of the file 

(you can type the first two to three letters of the file name, then just press "tab" and the file name should pop-up)

|
10.) The result files should be where-ever the file you ran was located.


