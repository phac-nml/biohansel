Degenerate Base Expansion
=========================

Biohansel allows users to input custom schemes containing `IUPAC degenerate bases <https://www.bioinformatics.org/sms/iupac.html>`_ in their kmers.
These degenerate bases allow for increased flexibility for kmer matches but can come at the cost of computing power.

In this section, we will look at the how to use degenerate bases and where to be careful when dealing with degenerate bases in
your biohansel genotyping schemes.

Including Degenerate Bases in k-mers:
-------------------------------------

Any of the degenerate bases can be added anywhere in the genotyping schemes k-mers. Once this happens, when the code is run, the k-mers 
are expanded into all possible combinations of normal DNA bases (A,G,C,T). These degenerate bases can be anywhere in the k-mer so long as there is
still a SNP separating the pair, and, if they are not in the SNP location the degenerate base matches in both k-mers of the pair.

Degenerate Base as a SNP:
#########################

Example k-mer pair with a single degenerate base as the SNP separating the pair. Here the SNP is bolded to try to make it easier to see. Remember
that biohansel scheme k-mers can be any length so long as they are the same length throughout the scheme:

| >1231-2.2
| TA\ **D**\ CT
|
| >negative1231-2.2
| TA\ **C**\ CT

The k-mer with the degenerate base would expand into three separate k-mers for the same genotype at the same position. Once the code was ran, 
these three positive k-mers would only count as 1 k-mer for the results outputs.

Expanded k-mers where if any of the positive k-mers are found in this set, you have genotype 2.2:

| >1231-2.2
| TA\ **A**\ CT
|
| >1231-2.2
| TA\ **T**\ CT
|
| >1231-2.2
| TA\ **G**\ CT
|
| >negative1231-2.2
| TA\ **C**\ CT


Multiple Degenerate Bases Elsewhere in the K-mer:
#################################################

As stated, you can have any number of degenerate bases in the SNP scheme k-mers. Here is an example of a pair containing 
two degenerate bases not in the position of the SNP. Here the degenerate bases are bolded along with their expansions:

| >1231-2.2
| CT\ **R**\ ACT\ **W**
|
| >negative1231-2.2
| CT\ **R**\ CCT\ **W**

This pair would expand into 4 different k-mers for each pair or **8** new k-mers in total. For the positive ones:

| >1231-2.2
| CT\ **A**\ ACT\ **A**
|
| >1231-2.2
| CT\ **A**\ ACT\ **T**
|
| >1231-2.2
| CT\ **G**\ ACT\ **A**
|
| >1231-2.2
| CT\ **G**\ ACT\ **T**

And for the negative ones:

| >negative1231-2.2
| CT\ **A**\ CCT\ **A**
|
| >negative1231-2.2
| CT\ **A**\ CCT\ **T**
|
| >negative1231-2.2
| CT\ **G**\ CCT\ **A**
|
| >negative1231-2.2
| CT\ **G**\ CCT\ **T**

As you can see, they are still only separated by the **A** SNP in the positive k-mer and the **C** SNP in the negative one.
Now there are 4 different pairs, with 8 separate sequences making up these 4 pairs. This expanded all from the
original pair and all of these new pairs would check for genotype 2.2 at position 1231.

In this example, there would end up being a total of 16 different k-mers due to the reverse compliment also being input into biohansel.


Unchecked Expansion of K-mers:
-------------------------------

When calculating the number of k-mers that are being created by a k-mer pair, the IUPAC DNA Nucleotide codes can be thought of
as the number of possibilities that they expand to and not what they represent as seen below:

| A, G, C, T = 1
|
| R, Y, S, W, K, M = 2
|
| B, D, H, V = 3
|
| N = 4

Using this, we can take any nucleotide k-mer we want and calculate the number k-mers that are being created.

Using the values above, we can substitute the nucleotide with a number and find the number of k-mers created.

Example sequence is "ACGTAGC":

| (A)(C)(G)(T)(A)(G)(C)
| (1)(1)(1)(1)(1)(1)(1) = 1

If we had a sequence with degenerate bases, we will start to see how fast the number of k-mers can increase.

Example sequence is "ACTNNANNTTA"

| (A) (C) (T) (N) (N) (A) (N) (N) (T) (T) (A)
| (1) (1) (1) (4) (4) (1) (4) (4) (1) (1) (1) = 256

This example is an example how having only 4 'N' nucleotides expands our one sequence into 256 different ones.

These 256 k-mers aren't including the second part of the pair that this sequence has to belong to and as their is a SNP of "A" in the middle, 
the other k-mer must also contain those 4 'N's. This means that there are 512 k-mers being used for just this pair alone.

Even then, the 512 k-mers for the positive and negative positions become a total of 1024 different k-mers due to the need to 
take into account the reverse compliment of all of the sequences.

To put this in perspective, the Heidelberg SNP genotyping scheme contains 202 pairs with 404 sequences and once ran, this is expanded
to 808 sequences by biohansel due to the reverse compliment input. The whole genotyping scheme has less k-mers than a single SNP pair in this case.

The goal is the remember that even a small number of degenerate bases can lead to a large number of k-mers and longer run times.
'N' is the extreme however and if you were creating a scheme with only the "2" value degenerate bases (ex. 'R'), then you could have
8 degenerate bases for a single pair and end up with the same 1024 expanded k-mers from the pair.

\* The take away here is to be careful when including degenerate bases in your scheme. The **more degenerate bases included**, the
**more kmers** are that are produced by expansion, the **slower the run time**, and the **more RAM is needed to run the sample**.


Benchmarking Degenerate Bases
-----------------------------

Expand Degenerate Base Module
#############################

More degenerate bases = More K-mers = Slower Run Times

In this section we are going to look at the speed of the expand base module and the code itself for different numbers of k-mers.

Here is the speed of running the **expand degenerate bases module** (not biohansel itself) on 1 core using pythons timeit:

+--------------------+-------------------------+----------------+  
| **K-mer Sequence** | **Max K-mers Produced** | **Time**       |
+--------------------+-------------------------+----------------+
| A                  | 1                       | 1.59 microsec  |
+--------------------+-------------------------+----------------+
| N                  | 4                       | 1.79 microsec  |
+--------------------+-------------------------+----------------+
| NN                 | 16                      | 2.47 microsec  |
+--------------------+-------------------------+----------------+
| NNN                | 64                      | 5.61 microsec  |
+--------------------+-------------------------+----------------+
| NNNN               | 256                     | 17.50 microsec |
+--------------------+-------------------------+----------------+
| NNNNN              | 1024                    | 68.30 microsec |
+--------------------+-------------------------+----------------+
| NNNNNN             | 4096                    | 305.0 microsec |
+--------------------+-------------------------+----------------+
| NNNNNNN            | 16384                   | 1.41 msec      |
+--------------------+-------------------------+----------------+
| NNNNNNNN           | 65536                   | 6.15 msec      |
+--------------------+-------------------------+----------------+
| NNNNNNNNN          | 262144                  | 26.5 msec      |
+--------------------+-------------------------+----------------+
| NNNNNNNNNN         | 1048576                 | 112.0 msec     |
+--------------------+-------------------------+----------------+
| NNNNNNNNNNN        | 4194394                 | 470.0 msec     |
+--------------------+-------------------------+----------------+
| NNNNNNNNNNNN       | 16777216                | 1.950 sec      |
+--------------------+-------------------------+----------------+
| NNNNNNNNNNNNN      | 67108864                | 8.930 sec      |
+--------------------+-------------------------+----------------+
| NNNNNNNNNNNNNN     | 268435456               | **Died**       |
+--------------------+-------------------------+----------------+

The higher k-mers are a bit of a stretch but show how much longer the module takes **PER K-MER** if you are not careful.

Remember, the above chart is for a **singular** k-mer and does not take into account the expansion of the whole scheme.
If you had a scheme with a lot of these, it would take that amount of time for each k-mer!


Whole Biohansel Code
####################

Benchmarking all of the biohansel code using the same input fasta file but increasing the total k-mer count each time.

Remember that the total number of k-mers if there are no degenerate bases is equal to the number of pairs multiplied by 
4 to take into account two sequences per pair and the RC of each sequence.

If degenerate bases are present, it is harder to guess the number and running biohansel will tell you if you have over the default 100,000 k-mers and
allow you to set the value that you deem acceptable with the "--max-degenerate-kmers" command.

+---------------------------+--------------------------------+--------------------------+  
| **Number of Nucleotides** | **Number of Scheme K-mers**    | **Execution Time (sec)** |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 808 --> Base Heidelberg Scheme | 0.613                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 4,394                          | 0.663                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 12,074                         | 0.721                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 36,650                         | 0.873                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 69,418                         | 0.971                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 134,954                        | 1.031                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 266,026                        | 1.502                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 528,171                        | 2.150                    |
+---------------------------+--------------------------------+--------------------------+
| 4751529                   | 1,052,459                      | 3.269                    |
+---------------------------+--------------------------------+--------------------------+

This work was done on an assembled fasta file. Not that even with 1,000,000 k-mers, the time it takes to run biohansel is only 3 seconds.
BUT, if you're using fastq files it is going to be much longer and they haven't been tested for speed with expansion yet! So be careful
with large expansions on fastq files.
