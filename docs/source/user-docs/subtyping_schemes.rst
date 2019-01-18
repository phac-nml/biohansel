
Subtyping Schemes 
================= 

.. |command| image:: https://raw.githubusercontent.com/phac-nml/biohansel/readthedocs/docs/source/user-docs/Screen%20Shot%202018-10-18%20at%203.22.52%20PM.png
   :alt: command line commands
   :width: 600 px

.. |scheme| image:: Genotype_scheme.png
   :alt: Genotype scheme chart
   :width: 600 px

This section will cover the subtyping schemes currently used by BioHansel for *Salmonella enterica* subspecies enterica serovar heidelberg and serovar enteritidis along with providing indepth information on how to create a custom subtyping scheme. 

Heidelberg and Enteritidis Subtyping Schemes 
--------------------------------------------  

The heidelberg subtyping scheme currently included with BioHansel features a set of 202 33-mer pairs with a single nucleotide polymorphism (SNP) distinguishing between the positive and negative condition allowing for the identification and classification of different subtypes based on the number and location of SNPs in a WGS sample. The enteritidis scheme features a set of 224 33-mer pairs in the same style as the heidelberg scheme to classify and identify different enteritidis serovar subtypes. Both schemes are maintained and developed by Genevieve Labbe et al. with changes occuring as new classifications are made.

The subtyping scheme created for each serovar follows a nested hierarchical approach to allow relationships between subtypes to be established based on which SNPs they contain. The SNP variations can link outbreaks and determine potential places for intervention. The format created additionally can be readily modified to take into account the constant slow evolution of the pathogen and provide places where new subtypes can be fit into the existing scheme.

The scheme and process that is used to subtype *Salmonella enterica* subspecies heidelberg can be seen below:

|scheme|

Following the scheme, each positive SNP match leads to a more specific classification of the subtype. All of the positive k-mers on the path to the final subtype should match the sample as seen in the diagram to get a positive subtype ID. If this does not occur, then BioHansel will output an error in the results. The `Output section <output.html>`_ contains more details on the errors that can occur.

The 33-mers that are found in the heidelberg.fasta file and the enteritidis.fasta included in the installation of BioHansel are based off of statistically conserved significant regions (no mobile elements) where the rate of SNP changes per year is small enough to isolate and relate subtypes.

The k-mers are structured as such:




