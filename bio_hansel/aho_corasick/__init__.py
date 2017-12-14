# -*- coding: utf-8 -*-

from collections import defaultdict

from ahocorasick import Automaton
import pandas as pd

from ..parsers import parse_fasta, parse_fastq
from ..utils import revcomp


def init_automaton(scheme_fasta):
    """Initialize Aho-Corasick Automaton with kmers from SNV scheme fasta

    Args:
        scheme_fasta: SNV scheme fasta file path

    Returns:
         Aho-Corasick Automaton with kmers loaded
    """
    A = Automaton()
    for h, s in parse_fasta(scheme_fasta):
        A.add_word(s, (h, s, False))
        A.add_word(revcomp(s), (h, s, True))
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
    for h, s in parse_fasta(fasta):
        for idx, (tilename, seq, is_revcomp) in A.iter(s):
            res.append((tilename, seq, is_revcomp, h, idx))
    df = pd.DataFrame(res, columns=['tilename', 'seq', 'is_revcomp', 'contig_id', 'match_index'])
    return df


def find_in_fastqs(A: Automaton, *fastqs):
    """Find scheme kmers in input fastq files

    Args:
        A: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fastqs: Input fastq file paths

    Returns:
        Dataframe with any matches found in input fastq files
    """
    res = defaultdict(int)
    for fq in fastqs:
        for h, s in parse_fastq(fq):
            for idx, (tilename, seq, is_revcomp) in A.iter(s):
                res[seq] += 1
    tmp = []
    for seq, freq in res.items():
        tilename, s, _ = A.get(seq)
        tmp.append((tilename, seq, freq))
    df = pd.DataFrame(tmp, columns=['tilename', 'seq', 'freq'])
    return df
