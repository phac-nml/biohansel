# -*- coding: utf-8 -*-
import logging
import re
from collections import defaultdict
from typing import Dict, List, Set

import attr

from .parsers import parse_fasta


@attr.s
class SubtypeCounts:
    subtype = attr.ib()
    refpositions = attr.ib(default=None)
    subtype_kmer_count = attr.ib(default=0, validator=attr.validators.instance_of(int))
    positive_kmer_count = attr.ib(default=0, validator=attr.validators.instance_of(int))
    negative_kmer_count = attr.ib(default=0, validator=attr.validators.instance_of(int))
    all_kmer_count = attr.ib(default=0, validator=attr.validators.instance_of(int))

    @subtype.validator
    def _check_subtype(self, attribute, value):
        REGEX_SUBTYPE = re.compile(r'^\d+(\.\d+)*$')
        if value is None or value == '':
            raise ValueError('Subtype cannot be None or empty string')
        if not REGEX_SUBTYPE.match(value):
            raise ValueError(f'Invalid subtype specified! The "{value}" kmer '
                             f'is not formatted correctly. It must be numbers '
                             f'delimited by "." (periods)')
        return value

    @subtype_kmer_count.validator
    def _check_subtype_kmer_count(self, attribute, value):
        if value == 0:
            raise ValueError('Subtype kmer count cannot be zero!')

    @positive_kmer_count.validator
    def _check_positive_kmer_count(self, attribute, value):
        if value < self.subtype_kmer_count:
            raise ValueError(f'Subtype {self.subtype}: Number of all subtype '
                             f'positive kmers (n={value}) cannot be less than the '
                             f'number of subtype specific kmers (n={self.subtype_kmer_count})')
        if value > self.all_kmer_count:
            raise ValueError(f'Subtype {self.subtype}: Number of all subtype '
                             f'positive kmers (n={value}) cannot exceed '
                             f'number of all matching kmers for subtype (n={self.subtype_kmer_count})')
        if value == 0:
            raise ValueError(f'Subtype {self.subtype}: Number of all subtype '
                             f'positive kmers cannot be zero!')


def _kmers(kmers_fasta: str) -> (Dict[str, List[str]], Dict[str, List[str]], Set[int]):
    kmers = defaultdict(list)
    neg_kmers = defaultdict(list)
    sizes = set()
    for h, s in parse_fasta(kmers_fasta):
        sizes.add(len(s))
        _, st = h.split('-')
        if 'negative' not in h:
            kmers[st].append(h)
        else:
            neg_kmers[st].append(h)
    return kmers, neg_kmers, sizes


def subtype_counts(scheme_fasta: str) -> Dict[str, SubtypeCounts]:
    subtype_counts = {}
    kmers, neg_kmers, sizes = _kmers(scheme_fasta)
    if len(sizes) > 1:
        logging.warning('Not all markers in "%s" of the same size! %s', scheme_fasta, sizes)
    ks = [x for x in kmers.keys()]
    ks.sort()
    for k in ks:
        st_pos_count_rest = 0
        subtypes_set = {k}
        if len(k) > 1:
            *r, _ = k.split('.')
            for i in range(len(r)):
                sub_k = '.'.join(r[:i + 1])
                subtypes_set.add(sub_k)
                pos_kmers = kmers[sub_k]
                st_pos_count_rest += len(pos_kmers)
        st_count = len(kmers[k])
        st_count_pos = st_pos_count_rest + st_count
        st_neg_count = sum(len(v) for k, v in neg_kmers.items() if k not in subtypes_set)

        subtype_count = SubtypeCounts(subtype=k,
                                      refpositions={int(v.split('-')[0]) for v in kmers[k]},
                                      subtype_kmer_count=st_count,
                                      positive_kmer_count=st_count_pos,
                                      negative_kmer_count=st_neg_count,
                                      all_kmer_count=st_neg_count + st_count_pos)
        subtype_counts[k] = subtype_count
    return subtype_counts
