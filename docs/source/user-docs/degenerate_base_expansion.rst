Degenerate Base Expansion
=========================

Biohansel allows users to input custom schemes containing `IUPAC degenerate bases <https://www.bioinformatics.org/sms/iupac.html>`_ in their kmers.
These degenerate bases allow for increased flexibility for kmer matches but can come at the cost of computing power.

In this section, we will look at the how to use degenerate bases and where to be careful when dealing with degenerate bases in
your biohansel subtyping schemes.

Including Degenerate Bases in k-mers:
-------------------------------------

Any of the degenerate bases can be added anywhere in the subtyping schemes k-mers. Once this happens, when the code is run, the k-mers 
are expanded into all possible combinations bases on the IUPAC bases. These degenerate bases can be anywhere in the k-mer so long as there is
still a SNP separating the pair, and, if they are not in the SNP location the degenerate base matches in both k-mers of the pair.

Degenerate Base as a SNP:
#########################

Example k-mer pair with a single degenerate base as the SNP separating the pair. Here the SNP is bolded to try to make it easier to see. Remember
that biohansel scheme k-mers can be any length so long as they are the same throughout the scheme:

| >1231-2.2
| TA\ **D**\ CT
|
| >negative1231-2.2
| TA\ **C**\ CT

The k-mer with the degenerate base would expand into three separate k-mers for the same subtype at the same position such 
that when ran, they only count as 1 k-mer and not three separate ones.

After base expansion, you would end up with three positive k-mers and one negative one. If any of the three positive ones are found,
then you have that subtype just as normal.

Expanded k-mers where if any of the positive k-mers are found you have subtype 2.2:

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

Example of a pair containing two degenerate bases not as a part of the SNP. Here the SNP is not bolded but the degenerate bases are and there
expansions will be too:

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
But now there are 4 different pairs, with 8 separate sequences that expanded from the original pair. With this, it shows just how fast the
k-mers number of k-mers can expand with only two degenerate bases. 

In this example, there would end up being a total of 16 different k-mers due to the reverse compliment also being input into biohansel.


Expansion of K-mer Numbers:
---------------------------

As demonstrated, the number of k-mers can grow rapidly and then is multiplied by two to take into account the reverse
compliment of the 