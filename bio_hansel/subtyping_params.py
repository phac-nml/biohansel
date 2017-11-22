# -*- coding: utf-8 -*-
import attr

@attr.s
class SubtypingParams(object):
    low_coverage_depth_freq = attr.ib(default=20, validator=attr.validators.instance_of(int))
    missing_total_tiles_max = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
    inconsistent_tiles_max = attr.ib(default=3, validator=attr.validators.instance_of(int))
    intermediate_subtype_tiles_max = attr.ib(default=0.05, validator=attr.validators.instance_of(float))
