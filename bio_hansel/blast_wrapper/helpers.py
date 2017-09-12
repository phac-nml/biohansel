# -*- coding: utf-8 -*-
import logging
import numpy as np

#: set: valid IUPAC nucleotide characters for checking FASTA format
VALID_NUCLEOTIDES = {'A', 'a',
                     'C', 'c',
                     'G', 'g',
                     'T', 't',
                     'R', 'r',
                     'Y', 'y',
                     'S', 's',
                     'W', 'w',
                     'K', 'k',
                     'M', 'm',
                     'B', 'b',
                     'D', 'd',
                     'H', 'h',
                     'V', 'v',
                     'N', 'n',
                     'X', 'x', }  # X for masked nucleotides

NT_SUB = {x: y for x, y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}


def revcomp(s):
    """Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return ''.join([NT_SUB[c] for c in s[::-1]])


def extend_subj_match_vec(df):
    """
    Get the extended clipped (clamped) start and end subject sequence indices

    Also get whether each match needs to be reverse complemented and whether each extended match would be truncated by
    the end of the subject sequence.

    Args:
        df (pandas.DataFrame): blastn results dataframe

    Returns:
        int pandas.Series: extended and clipped start indices
        int pandas.Series: extended and clipped end indices
        bool pandas.Series: does extracted seq need to be reverse complemented?
        bool pandas.Series: would the extended seq be truncated by the ends of the subject sequence?
        bool pandas.Series: was the subject seq extended?
    """
    needs_revcomp = df.sstart > df.send
    add_to_end = df.qlen - df.qend
    add_to_start = df.qstart - 1
    ssum2 = (df.send + df.sstart) / 2.0
    sabs2 = np.abs(df.send - df.sstart) / 2.0
    end_idx = ssum2 + sabs2 - 1
    start_idx = ssum2 - sabs2 - 1
    start_idx[needs_revcomp] -= add_to_end
    start_idx[~needs_revcomp] -= add_to_start
    end_idx[needs_revcomp] += add_to_start
    end_idx[~needs_revcomp] += add_to_end
    clipped_start_idx = np.clip(start_idx, 0, (df.slen - 1))
    clipped_end_idx = np.clip(end_idx, 0, (df.slen - 1))
    trunc = (clipped_start_idx != start_idx) | (clipped_end_idx != end_idx)
    is_extended = (add_to_start > 0) | (add_to_end > 0)
    return clipped_start_idx, clipped_end_idx, needs_revcomp, trunc, is_extended


def retrieve_seq(seq, start, end, needs_revcomp):
    """Retrieve sub-sequence from nucleotide sequence

    Args:
        seq (str): nucleotide sequence
        start (int): Start index
        end (end): End index
        needs_revcomp (bool): reverse complement subseq?

    Returns:
         str: subseq of `seq`
    """
    start = int(start)
    end = int(end)
    out_seq = seq[start:(end + 1)]
    if needs_revcomp:
        out_seq = revcomp(out_seq)
    return out_seq


def parse_fasta(filepath):
    '''Parse a fasta file returning a generator yielding tuples of fasta headers to sequences.

    Note:
        This function should give equivalent results to SeqIO from BioPython

        .. code-block:: python

            from Bio import SeqIO
            # biopython to dict of header-seq
            hseqs_bio = {r.description:str(r.seq) for r in SeqIO.parse(fasta_path, 'fasta')}
            # this func to dict of header-seq
            hseqs = {header:seq for header, seq in parse_fasta(fasta_path)}
            # both methods should return the same dict
            assert hseqs == hseqs_bio

    Args:
        filepath (str): Fasta file path

    Returns:
        generator: yields tuples of (<fasta header>, <fasta sequence>)
    '''
    with open(filepath, 'r') as f:
        seqs = []
        header = ''
        line_count = 0
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line[0] == '>':
                if header == '':
                    header = line.replace('>', '')
                else:
                    yield header, ''.join(seqs)
                    seqs = []
                    header = line.replace('>', '')
            else:
                non_nucleotide_chars_in_line = set(line) - VALID_NUCLEOTIDES
                if len(non_nucleotide_chars_in_line) > 0:
                    msg = '{file}: Line {line} contains the following non-nucleotide characters: {chars}'.format(
                        file=filepath,
                        line=line_count,
                        chars=', '.join([x for x in non_nucleotide_chars_in_line]))
                    logging.warning(msg)
                seqs.append(line)
            line_count += 1
        yield header, ''.join(seqs)
