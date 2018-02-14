# -*- coding: utf-8 -*-
import attr

@attr.s
class SubtypingParams(object):
    low_coverage_depth_freq = attr.ib(default=20.0, validator=attr.validators.instance_of((float, int)))
    max_perc_missing_tiles = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    min_ambiguous_tiles = attr.ib(default=3, validator=attr.validators.instance_of(int))
    max_perc_intermediate_tiles = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    min_kmer_freq = attr.ib(default=8, validator=attr.validators.instance_of((float, int)))
    max_kmer_freq = attr.ib(default=1000, validator=attr.validators.instance_of((float, int)))
    min_coverage_warning = attr.ib(default=20, validator=attr.validators.instance_of((float, int)))

    @max_perc_missing_tiles.validator
    def _validate_max_perc_missing_tiles(self, attribute, value):
        if 0.0 > value > 1.0:
            raise AttributeError('Max % misssing tiles was {} expected to be decimal between 0.0 and 1.0 inclusive'.format(value))
