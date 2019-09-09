# -*- coding: utf-8 -*-

from collections import defaultdict
from itertools import product
import logging
import sys

from ahocorasick import Automaton
import pandas as pd

from ..parsers import parse_fasta, parse_fastq
from ..utils import revcomp
from ..const import bases_dict
from ..subtyping_params import SubtypingParams

def expand_degenerate_bases(seq):
    """List all possible kmers for a scheme given a degenerate base
   
    Args:
         Scheme_kmers from SNV scheme fasta file

    Returns:
         List of all possible kmers given a degenerate base or not
    """

    return list(map("".join, product(*map(bases_dict.get, seq))))

def check_total_kmers(scheme_fasta, subtyping_params):
    '''Checks that the number of kmers about to be created is not at too high a computation or time cost

    Args:
         scheme_fasta: Kmer sequences from the SNV scheme
         Subtyping parameter max_degenerate_kmers: The max kmers allowed by the scheme
    
    Returns:
         None if created kmers < max degenerate kmers argument
         Exits code with warning if created kmers > max degenerate kmers argument
    '''
    kmer_number = 0
    for header, sequence in parse_fasta(scheme_fasta):
        value = 1
        for char in sequence:
            length_key = len(bases_dict[char])
            value = value * length_key
        kmer_number = kmer_number + value
    if kmer_number*2 > subtyping_params.max_degenerate_kmers:
        return logging.error(
            '''
    Your current scheme contains "{}" kmers which is over the current max-degenerate-kmers check of "{}" (Maximum recommended k-mers is 100000).
    It is not advised to run this scheme due to the time and memory usage required to give an output with this many kmers loaded.
    If you still want to run this scheme, add the command line check of "--max-degenerate-kmers {}" at the end of your previous command.
            '''.format(
                kmer_number*2,
                subtyping_params.max_degenerate_kmers,
                kmer_number*2+1
                )), sys.exit()
    else:
        return None


def init_automaton(scheme_fasta):
    """Initialize Aho-Corasick Automaton with kmers from SNV scheme fasta

    Args:
        scheme_fasta: SNV scheme fasta file path

    Returns:
         Aho-Corasick Automaton with kmers loaded
    """
    A = Automaton()
    for header, sequence in parse_fasta(scheme_fasta):
        kmer_list = expand_degenerate_bases(sequence)
        for seq in kmer_list:
            A.add_word(seq, (header, seq, False))
            A.add_word(revcomp(seq), (header, seq, True))
    A.make_automaton()
    return A


def find_in_fasta(A: Automaton, fasta: str) -> pd.DataFrame:
    """Find scheme kmers in input fasta file

    Args:
        A: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fasta: Input fasta path

    Returns:
        Dataframe with any matches found in input fasta file
    """
    res = []
    for contig_header, sequence in parse_fasta(fasta):
        for idx, (kmername, kmer_seq, is_revcomp) in A.iter(sequence):
            res.append((kmername, kmer_seq, is_revcomp, contig_header, idx))
    df = pd.DataFrame(res, columns=['kmername', 'seq', 'is_revcomp', 'contig_id', 'match_index'])
    return df


def find_in_fastqs(A: Automaton, *fastqs):
    """Find scheme kmers in input fastq files

    Args:
        A: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fastqs: Input fastq file paths

    Returns:
        Dataframe with any matches found in input fastq files
    """
    kmer_seq_counts = defaultdict(int)
    for fastq in fastqs:
        for _, sequence in parse_fastq(fastq):
            for idx, (_, kmer_seq, _) in A.iter(sequence):
                kmer_seq_counts[kmer_seq] += 1
    res = []
    for kmer_seq, freq in kmer_seq_counts.items():
        kmername, sequence, _ = A.get(kmer_seq)
        res.append((kmername, kmer_seq, freq))
    df = pd.DataFrame(res, columns=['kmername', 'seq', 'freq'])
    return df
