# -*- coding: utf-8 -*-
import attr


@attr.s
class SubtypingParams(object):
    low_coverage_depth_freq = attr.ib(default=20.0, validator=attr.validators.instance_of((float, int)))
    max_perc_missing_kmers = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    min_ambiguous_kmers = attr.ib(default=3, validator=attr.validators.instance_of(int))
    max_perc_intermediate_kmers = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    min_kmer_freq = attr.ib(default=8, validator=attr.validators.instance_of((float, int)))
    min_kmer_frac = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    max_kmer_freq = attr.ib(default=1000000, validator=attr.validators.instance_of((float, int)))
    min_coverage_warning = attr.ib(default=20, validator=attr.validators.instance_of((float, int)))
    max_degenerate_kmers = attr.ib(default=10000000, validator=attr.validators.instance_of(int))

    @max_perc_missing_kmers.validator
    def _validate_max_perc_missing_kmers(self, attribute, value):
        if 0.0 > value > 1.0:
            raise AttributeError(f'Max % misssing kmers was {value} expected '
                                 f'to be decimal between 0.0 and 1.0 inclusive')
