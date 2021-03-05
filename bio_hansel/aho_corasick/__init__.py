# -*- coding: utf-8 -*-

from collections import defaultdict

import pandas as pd
from ahocorasick import Automaton

from ..parsers import parse_fasta, parse_fastq
from ..utils import revcomp, expand_degenerate_bases


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


def find_in_fasta(automaton: Automaton, fasta: str) -> pd.DataFrame:
    """Find scheme kmers in input fasta file

    Args:
        automaton: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fasta: Input fasta path

    Returns:
        Dataframe with any matches found in input fasta file
    """
    res = []
    for contig_header, sequence in parse_fasta(fasta):
        for idx, (kmername, kmer_seq, is_revcomp) in automaton.iter(sequence):
            res.append((kmername, kmer_seq, is_revcomp, contig_header, idx))
    columns = ['kmername', 'seq', 'is_revcomp', 'contig_id', 'match_index']
    return pd.DataFrame(res, columns=columns)


def find_in_fastqs(automaton: Automaton, *fastqs):
    """Find scheme kmers in input fastq files

    Args:
        automaton: Aho-Corasick Automaton with scheme SNV target kmers loaded
        fastqs: Input fastq file paths

    Returns:
        Dataframe with any matches found in input fastq files
    """
    kmer_seq_counts = defaultdict(int)
    for fastq in fastqs:
        for _, sequence in parse_fastq(fastq):
            for idx, (_, kmer_seq, _) in automaton.iter(sequence):
                kmer_seq_counts[kmer_seq] += 1
    res = []
    for kmer_seq, freq in kmer_seq_counts.items():
        kmername, sequence, _ = automaton.get(kmer_seq)
        res.append((kmername, kmer_seq, freq))
    return pd.DataFrame(res, columns=['kmername', 'seq', 'freq'])
