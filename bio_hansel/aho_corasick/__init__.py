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
    for header, sequence in parse_fasta(scheme_fasta):
        A.add_word(sequence, (header, sequence, False))
        A.add_word(revcomp(sequence), (header, sequence, True))
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
        for idx, (tilename, tile_seq, is_revcomp) in A.iter(sequence):
            res.append((tilename, tile_seq, is_revcomp, contig_header, idx))
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
    tile_seq_counts = defaultdict(int)
    for fastq in fastqs:
        for _, sequence in parse_fastq(fastq):
            for idx, (_, tile_seq, _) in A.iter(sequence):
                tile_seq_counts[tile_seq] += 1
    res = []
    for tile_seq, freq in tile_seq_counts.items():
        tilename, sequence, _ = A.get(tile_seq)
        res.append((tilename, tile_seq, freq))
    df = pd.DataFrame(res, columns=['tilename', 'seq', 'freq'])
    return df
