# -*- coding: utf-8 -*-
import attr

@attr.s
class Subtype(object):
    sample = attr.ib(validator=attr.validators.instance_of(str))
    file_path = attr.ib(validator=attr.validators.instance_of(str))
    scheme = attr.ib(validator=attr.validators.instance_of(str))
    scheme_version = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    subtype = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    all_subtypes = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    inconsistent_subtypes = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    tiles_matching_subtype = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    are_subtypes_consistent = attr.ib(default=True, validator=attr.validators.instance_of(bool))
    n_tiles_matching_all = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_positive = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_subtype = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_all_expected = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_positive_expected = attr.ib(default=0, validator=attr.validators.instance_of(int))
    n_tiles_matching_subtype_expected = attr.ib(default=0, validator=attr.validators.instance_of(int))
