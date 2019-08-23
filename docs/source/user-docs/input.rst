Input
=====

This section describes the main input files that are needed to run biohansel along with the various
methods available to analyze datasets using the tool.

The three input files are:

- The WGS data/assembled genome added as a FASTQ or FASTA file

- A chosen genotyping scheme (heidelberg or enteritidis) or a user created custom genotyping scheme (FASTA).
More info in the `genotyping schemes section <genotyping_schemes.html>`_.

- A metadata table (**Optional**) in CSV or TSV (tab-delimited) format to add additional information to the results. A **.tsv file is highly recommended** as csv files are currently unstable

More detailed information on the output of what each results file contains can be found in the `Output section <output.html>`_.


Types of Analysis
#################

Analysis of a single FASTA file
-------------------------------

Analysis of a single FASTA file would be run on a sample that has already been assembled into contigs
through another program/tool. To run the command, the arguments that must be specified include:

- -s "scheme" where the scheme defined can be one of the five built in schemes or a user created one (FASTA format)

- Any combination of the results delimiters and names (file names can be changed but must be included after the argument):

    - -o results.tab
    - -O match_results.tab
    - -S tech_results.tab

- The name of the FASTA file at the end of the command

An example command for the analysis of a single FASTA file called SRR1002850.fasta would look like:

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab /path/to/SRR1002850.fasta

Or, if you have already changed to the directory containing the dataset, you can use the following
command where you do not have to specify the path to the data:

.. code-block:: bash

    hansel -s heidelberg -vv -o results.tab -O match_results.tab SRR1002850.fasta

The output of the biohansel tool can be found in the directory that the command was run from.



Analysis of a single FASTQ readset
----------------------------------

Analysis of a single FASTQ readset would be run on raw sequencing data. To run the command, the arguments that must be specified include:

- -s "scheme" where the scheme defined can be one of the five built in schemes or a user created one (FASTA format)

- Any combination of the results delimiters and names (file names can be changed but must be included after the argument):
 
    - -o results.tab
    - -O match_results.tab
    - -S tech_results.tab

- The name of the FASTQ file(s) at the end of the command
    - For single-end reads include the one file 
    - For paired-end reads include: -p followed by both files one after the other

An example command for the analysis of a single single-end reads run dataset would look like:

.. code-block:: bash

    hansel -s heidelberg -vv -t 4 -o results.tab -O match_results.tab SRR5646583.fastqsanger


An example command for the analysis of a single paired-end reads run dataset would look like:

.. code-block:: bash

    hansel -s heidelberg -vv -t 4 -o results.tab -O match_results.tab -p SRR5646583_forward.fastqsanger SRR5646583_reverse.fastqsanger


Analysis of all FASTA/FASTQ files in a directory
------------------------------------------------

Analysis on **all** of the FASTA/FASTQ files in the specified directory. This will run on all FASTA/FASTQ files
in the directory. Be sure that there are no miscellaneous files that may unnecessarily increase analysis time or
lead to unneeded errors.

Biohansel will only attempt to analyze the FASTA/FASTQ files within the specified directory and will
not descend into any subdirectories! As such, make sure all of the data to be analyzed is in the same
location or organized in a way that suits the project.

Analysis of all of the sequencing files in a directory must include following the arguments to run properly:

- -s "scheme" where the scheme defined can be one of the five built in schemes or a user created one (FASTA format)

- Any combination of the results delimiters and names (file names can be changed but must be included after the argument):
 
    - -o results.tab
    - -O match_results.tab
    - -S tech_results.tab

- -D /path/to/directory_with_data

Optionally, you are able to specify the number of threads for an analysis with the
--threads argument. If you do not specify this, it will default to 1.

- --threads <#_cpu> to specify the number of CPUs wanted to run the analysis.

An example of a general command for the analysis of a directory of FASTA/FASTQ files:

.. code-block:: bash

    hansel -s heidelberg -vv --threads <n_cpu> -o results.tab -O match_results.tab -D /path/to/fastas_or_fastqs/

The chosen output files can be found in the directory that the command was run from or that was specified in the output
names and it will contain data from each of the analyzed files run by biohansel. 

Ex. If I was running an analysis on samples stored in my "data" directory found in the path science/user/data,
I could cd to my user folder and run the following command:

.. code-block:: bash

    hansel -s heidelberg -vv --threads 1 -o results.tab -O match_results.tab -D data/


Genotype Metadata Table (Optional)
##################################

Optionally you can select a genotype metadata information table to include genotype metadata along with the genotyping
results created with biohansel. Metadata tables must be in a tab-delimited format to correctly work. The file extension
for your metadata table should be **.tsv** if at all possible or you may end up with an error and no analysis results.

To add a metadata table to the analysis you will add the argument `-M <metadata_scheme.tsv>` to any other analysis command.
There are no requirements for the number of columns or the content of each of the columns on the metadata table so
long as the **first column** is labeled as **"subtype"**. 

A command that incorporates the -M command for analysis would be structured following the previously established requirements and looks as follows:

.. code-block:: bash

    hansel -s heidelberg -M <metadata_scheme.tsv> -vv -o results.tab -O match_results.tab <data>


The biohansel results table will be joined with the genotype metadata table based if a genotype on the metadata
table matches one on the results. If a match occurs, the metadata of that genotype will be added to the table
at the end of the results.tab and tech_results.tab results files. 

Example metadata table (called metadata.tsv):

+-------------+-------+--------+----------+ 
| subtype     | Clade | Source | Symptoms | 
+=============+=======+========+==========+  
| 1           | I     | Geese  | Death    | 
+-------------+-------+--------+----------+ 
| 1.1         | I     | Moose  | Burns    | 
+-------------+-------+--------+----------+ 
| 2.2.1.1.1   | II    | Mouse  | Boils    | 
+-------------+-------+--------+----------+  
| 2.2.2.2.2.1 | IIa   | Human  | Rash     | 
+-------------+-------+--------+----------+ 

***When naming a metadata table make sure there are no spaces or parentheses and that its extension is .tsv or the analysis may fail.*** 

The added metadata will appear at then end of the results.tab and the tech_results.tab files.

Example: tech_results.tab without metadata added:

+----------+-----------+-----------+------------+ 
| sample   | subtype   | qc_status | qc_message | 
+==========+===========+===========+============+  
| CP012921 | 2.2.3.1.2 | PASS      |            | 
+----------+-----------+-----------+------------+ 


Example: tech_results.tab with metadata:

+----------+-----------+-----------+------------+-------+--------+----------+ 
| sample   | subtype   | qc_status | qc_message | Clade | Source | Symptoms | 
+==========+===========+===========+============+=======+========+==========+  
| CP012921 | 2.2.3.1.2 | PASS      |            | I     | Geese  | Rash     | 
+----------+-----------+-----------+------------+-------+--------+----------+ 


You can add metadata to the analysis with Galaxy by uploading either a .tsv or a .csv
file to your history and specifying that you want it used in the analysis. **A .tsv file is recommended**.
