
Subtyping Schemes 
================= 

.. |scheme| image:: Genotype_scheme.png
   :alt: Genotype scheme chart
   :width: 600 px

This section will cover the subtyping schemes currently used by BioHansel for *Salmonella enterica* subspecies enterica serovar heidelberg and serovar enteritidis along with providing in depth information on how to create a custom subtyping scheme. 

Heidelberg and Enteritidis Subtyping Schemes 
--------------------------------------------  

The heidelberg subtyping scheme currently included with BioHansel features a set of 202 33-mer pairs with a single nucleotide polymorphism (SNP) distinguishing between the positive and negative condition in the pair. This distinction between pairs allows for the identification and classification of different subtypes of heidelberg serovars based on the number and location of SNPs in a WGS sample. The enteritidis scheme features a similar set of 224 33-mer pairs in the same style as the heidelberg scheme to classify and identify different enteritidis serovar subtypes. Both schemes were developed by Genevieve Labbe et al. with changes occurring as new classifications are made.

Subtype Classification System
-----------------------------

The subtyping classification system created for BioHansel follows a nested hierarchical approach to allow relationships between subtypes to be established based on which SNPs they contain. The format designed for the classification system allows modification of the existing subtyping scheme to allow new branches to be defined as new subtypes are fit into the existing classification system. The also works as a way to easily link outbreak origins and look at places further up the hierarchy where interventions can be done. 

The scheme and process that is used to subtype *Salmonella enterica* subspecies heidelberg can be seen below:

|scheme|

Following the scheme, each positive k-mer match leads to a more specific classification of the subtype with the first matches determining which lineage the isolate is from (1 or 2 in this case). After the main lineage is determined, the subtype is determined based on all of the matching positive k-mers found in the sample as it follows along the path to the specific subtype. It is important that all/most positive and negative k-mers should match a spot in the sample to allow correct subtyping!

The `Output section <output.html>`_ contains more details on the errors that can be run into when running a sample.

Creating a Subtyping Scheme
---------------------------

Creating a representative and well established subtyping scheme for a BioHansel is a large task; but once it is established, a created scheme is easy to modify to fit the needs of the research and allow for new classifications as they are discovered. When creating a subtyping scheme, keep in mind that the organism should be highly clonal as all of the k-mers identified and created should be found in all/almost all isolates for BioHansel to work correctly. 

Detailed Steps
##############

The detailed steps to create a well structured and accurate subtyping scheme are as follows. These steps were used to create the Heidelberg and Enteritidis Subtyping Schemes and have been shown to create accurate results from the test samples run. The steps are:

1. Generate a large data set that is representative of the organisms population being defined. For best results make sure to:

- Remove outliers

- Remove poor quality data

- de-duplicate the data set

2. Choose an available reference genome for the organism (ideally closed)

3. Subdivide the population into closely related clonal groups using MASH followed by SNP analysis

- Aim for groups that are less then 3000SNPs between strains over more then 80% of the reference genome

4. Remove rare outliers from the data set

- these are detected by SNP matrices, number of unaligned bases, number of heterozygous sites, number of bases with low coverage, etc.

- These rare outliers are from suspected poor quality WGS data, mixed culture samples, or large recombinant regions (phage or transposons).

5. De-duplicate the data once again by removing strains that are nearly identical to each other. This can be defined as:

- Strains that are 0-2 SNPs apart over more then 80% of the reference genome

- Strains that MASH cluster with a distance of ≤ 0.001

6. Create a Maximum Likelihood (ML) phylogenetic tree from the SNP derived reference assembly of the strains to the reference genome. Here you are looking for:

- Regions that are conserved across the whole population of interest such that the SNPs in the areas are found in 99.5% of all isolates

- SNPs that are at least 20 base pairs from other SNPs or indels.

	- **The 20 bases on either side of the SNP should be conserved in at least 99.5% of isolates!**

7. Divide the ML tree into main lineages and sub-lineages according to the shape of the tree to allow users to identify the main clonal expansions. When doing this make sure that:

- Tree branches are at least 2 SNPs long

	- Longer the branch the better as there will be more SNP positions to choose from for defining that subtype

If wanted, you can lower the number of SNP sites to be evaluated into the scheme by removing all of the SNPs that are present in less then 5 isolates and then remake the tree. The aim is to have a least 5-10 strains per sub-lineage, to keep the scheme focused on clonal expansions.

8. Create a neighbour-joining tree and root it using a distantly related sequence or a pseudo sequence to determine where the root of the tree should be.

|
9. Give main lineages and sub-lineages determined previously hierarchical codes based on how they cluster in the NJ tree and the SNPs that make up each sequence.

|
10. Extract from the SNV table or VCF file the cannonical SNPs that define the subtype and differentiate it from other strains using `FEHT <https://github.com/chadlaing/feht>`_ which can be installed into bioconda or galaxy. 

The installation instructions are found in the link but if you are using bioconda for BioHansel, the easiest thing to do is go to the correct environment and install FEHT there with the following commands:

.. code-block:: bash

    conda activate <name of environment to install feht to>

    conda install -c bioconda feht

FEHT needs the following specific files to run this process:

- A metadata file with the hierarchical codes

- A SNV table or a VCF file that defines the subtype

**Make sure that the isolate names match exactly and both files use a tab delimiter**

The metadata file should look as such and be in a **.tsv** format:

+---------------+---------+---------+---------+---------+-----+
| Strain name   | Level 1 | Level 2 | Level 3 | Level 4 | ... |
+===============+=========+=========+=========+=========+=====+  
| SRR1242421444 | 1       | 1.1     | 1.1.2   | 1.1.2.3 | ... |
+---------------+---------+---------+---------+---------+-----+  
| SRR1242422313 | 2       | 2.2     | 2.2.2   | 2.2.2   | ... |
+---------------+---------+---------+---------+---------+-----+

11. Extract the exact matches to the query using the ratioFilter in FEHT by switching "-f" to "1". 

This is done as the FEHT program performs an all-against-all comparison of all the subtypes, one column (one hierarchy) at a time and we only want the exact matches.

|
12. From this output, we want to extract the subtype against all else results by searching for the ! sign (ex. search !2.2 instead of 2.2) and compile these results into a new **.tsv** file with the following information:

+---------+--------------+---------------+---------------+
| Subtype | SNP Location | Positive Base | Negative Base |
+=========+==============+===============+===============+
| 1       | 395          | A             | G             | 
+---------+--------------+---------------+---------------+
| 1       | 2998         | T             | G             | 
+---------+--------------+---------------+---------------+
| 1.1     | 29231        | A             | G             | 
+---------+--------------+---------------+---------------+
| 1.1.1   | 77889        | T             | C             | 
+---------+--------------+---------------+---------------+

The positive base is the base found in the middle of the k-mer and it corresponds to the subtype of the sample. The negative base is the base found in all other samples. Both are equally important for the program to function properly so it is essential that they are properly defined.

13. Create the subtyping scheme with all of the information obtained. The SNP column shows the exact position that the SNP is found in the reference genome. This spot can be made into a 33-mer tile used in the scheme by recording 16 bases on each side of the SNP such that the SNP is in position 17 of the 33-mer.

A script can be used to do this which will create 33-mers from the reference genome. Keep in mind that most of them will be of the negative variety and the positive k-mer pair will need to be created.

14. Finish the subtyping scheme by making sure that each carefully crafted 33-mer has a positive and negative pair attatched to the correct subtype. This can be done also using a script (currently being worked on) or the following method:


    1. Paste the 33-mers into the correct location in the FEHT filtered output spreadsheet next to the corresponding SNPs.  

    2. The 33 bp sequences are expanded using TextWrangler (replace [A,T,C,G] by the same base+tab), then pasted back into excel, in 33 adjacent columns.  

    3. Replace the 17th column (middle one) with the positive base column, and collapse the 33 columns into one by removing the tabs in text wrangler.  

    4. Paste back into Excel as the list of “positive tiles”.  

    5. Replace the middle column by the negative base column and repeat the same procedure to obtain the list of “negative tiles”.

15. Create a FASTA file following the K-mer structure found below. Make sure that the headers and sequences are on separate lines. The order of the files in the scheme does not matter for BioHansel input.

It is important that the K-mers follow the exact format or the analysis will generate errors and potentially fail. They should all be the same size with position 17 (or the middle position) containing the SNP.

K-mer Structure
###############

The structure k-mer pairs are structured as such and must follow the following format:

| >[SNP position in ref genome]-[subtype] for the positive tiles
| AAATTTCAGCTAGCTAGCTAGCAATCACTGATC
| 
| >negative[SNP position in ref genome]-[subtype] for the negative tiles
| AAATTTCAGCTAGCTATCTAGCAATCACTGATC

An example with real data:

| >2981-2.2.3.1.4
| ACTGCCGCCGGAGCCGTGTGAAAATATTGTTTA
| 
| >negative2981-2.2.3.1.4
| ACTGCCGCCGGAGCCGCGTGAAAATATTGTTTA


***The first distinction between subtypes 1 and 2 (or potentially more subtypes) does not have a negative condition and instead moves samples into one of the two classes established. The setup for the k-mers is similar to the other k-mers shown above:

| >717-1
| ATGCAGAGTCAGTCAGATCAACATGCACCCACA
| 
| >717-2
| ATGCAGAGTCAGTCAGTTCAACATGCACCCACA

16. Test the created scheme by running BioHansel to verify that all of the expected positive target sequences are present in the corresponding strains. Eliminate targeted k-mers from the scheme that do not work well and verify that the targeted k-mers created are present in most of the data set. Finally test the scheme on a de novo assembly along with raw Illumina sequencing reads to make sure it holds true for both.


