Degenerate Base Expansion
=========================

Biohansel allows users to input custom schemes containing `IUPAC degenerate bases <https://www.bioinformatics.org/sms/iupac.html>`_ in their kmers.
These degenerate bases allow for increased flexibility for kmer matches but can come at the cost of computing power.

In this section, we will look at the how to use degenerate bases and where to be careful when dealing with degenerate bases in
your biohansel subtyping schemes.

Including Degenerate Bases in k-mers:
-------------------------------------

Any of the degenerate bases can be added anywhere in the subtyping schemes k-mers. Once this happens, when the code is run, the k-mers 
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

The k-mer with the degenerate base would expand into three separate k-mers for the same subtype at the same position. Once the code was ran, 
these three positive k-mers would only count as 1 k-mer for the results outputs.

Expanded k-mers where if any of the positive k-mers are found in this set, you have subtype 2.2:

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
original pair and all of these new pairs would check for subtype 2.2 at position 1231.

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

This example is an example how having only 4 'N' nucleotides expands our one sequence into 256.

These 256 k-mers aren't including the second part of the pair that this sequence has to belong to and as their is a SNP A in the middle, 
the other k-mer must also contain those 4 'N's. This means that there are 512 k-mers being used for just this pair

The 512 k-mers become 1024 due to taking into account the reverse compliment of all of the sequences.

To put this in perspective, the Heidelberg SNP subtyping scheme contains 202 pairs with 404 sequences and once ran, this is expanded
to 808 sequences by biohansel due to reverse compliments. The whole scheme has less k-mers than a single SNP pair in this case.

The goal is the remember that even a small number of degenerate bases can lead to a large number of k-mers and longer run times.
'N' is the extreme however and if you were creating a scheme with only the "2" value degenerate bases (ex. 'R'), then you could have
8 degenerate bases and end up with the same 1024 expanded k-mers from the pair.

\* The take away is to be careful when including degenerate bases in your scheme. The **more degenerate bases included**, the **more kmers** are 
produced, the **slower the run time**, and the **more RAM is needed to run the sample**.


Benchmarking Degenerate Bases
-----------------------------

More degenerate bases = More K-mers = Slower Run Times

Even with the risks, they are included as they are useful. Here is the speed of running the **expand degenerate** 
**base module** (not biohansel itself) using pythons timeit module and different k-mers.


