# -*- coding: utf-8 -*-
import attr

@attr.s
class SubtypingParams(object):
    low_coverage_depth_freq = attr.ib(default=20, validator=attr.validators.instance_of(int))
    max_perc_missing_tiles = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    min_ambiguous_tiles = attr.ib(default=3, validator=attr.validators.instance_of(int))
    max_perc_intermediate_tiles = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    calc_min_kmer_freq = attr.ib(default=False, validator=attr.validators.instance_of(bool))
    min_kmer_freq = attr.ib(default=10, validator=attr.validators.instance_of(int))
    max_kmer_freq = attr.ib(default=200, validator=attr.validators.instance_of(int))
    savgol_window_len = attr.ib(default=3, validator=attr.validators.instance_of(int))
    savgol_poly_degree = attr.ib(default=2, validator=attr.validators.instance_of(int))
    kmer_cov_perc_confidence = attr.ib(default=0.95, validator=attr.validators.instance_of(float))

    @max_perc_missing_tiles.validator
    def _validate_max_perc_missing_tiles(self, attribute, value):
        if 0.0 > value > 1.0:
            raise AttributeError('Max % missing tiles was {} expected to be decimal between 0.0 and 1.0 inclusive'
                                 ''.format(value))

    @max_perc_intermediate_tiles.validator
    def _validate_max_perc_intermediate_tiles(self, attribute, value):
        if 0.0 > value > 1.0:
            raise AttributeError('Max % intermediate tiles was {} '
                                 'expected to be decimal between 0.0 and 1.0 inclusive'.format(value))

    @kmer_cov_perc_confidence.validator
    def _validate_kmer_cov_perc_confidence(self, attribute, value):
        if 0.0 > value > 1.0:
            raise AttributeError('Kmer Cov % Confidence was {} '
                                 'expected to be decimal between 0.0 and 1.0 inclusive'.format(value))
