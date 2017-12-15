# -*- coding: utf-8 -*-
from typing import List, Optional

import attr

from .const import REGEX_FASTQ


@attr.s
class Subtype(object):
    sample = attr.ib(validator=attr.validators.instance_of(str))
    file_path = attr.ib(validator=attr.validators.instance_of(str))
    scheme = attr.ib(validator=attr.validators.instance_of(str))
    scheme_version = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    subtype = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    non_present_subtypes = attr.ib(default=None)  # type: Optional[List[str]]
    all_subtypes = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    inconsistent_subtypes = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    tiles_matching_subtype = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    negative_tiles_matching_subtype = attr.ib(default=None,
                                              validator=attr.validators.optional(attr.validators.instance_of(str)))
    are_subtypes_consistent = attr.ib(default=True, validator=attr.validators.instance_of(bool))
    n_tiles_matching_all = attr.ib(default=0)
    n_tiles_matching_positive = attr.ib(default=0)
    n_tiles_matching_negative = attr.ib(default=0)
    n_tiles_matching_subtype = attr.ib(default=0)
    n_tiles_matching_all_expected = attr.ib(default=0)
    n_tiles_matching_positive_expected = attr.ib(default=0)
    n_tiles_matching_negative_expected = attr.ib(default=0)
    n_tiles_matching_subtype_expected = attr.ib(default=0)
    n_negative_tiles_matching_subtype_expected = attr.ib(default=0)
    qc_status = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    qc_message = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))
    scheme_subtype_counts = attr.ib(default=None, repr=False)
    input_type = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))

    def is_fastq_input(self):
        if isinstance(self.file_path, str):
            return bool(REGEX_FASTQ.match(self.file_path))
        elif isinstance(self.file_path, list):
            return all([bool(REGEX_FASTQ.match(x)) for x in self.file_path])
        else:
            raise ValueError('Unexpected type "{}" for "file_path": {}'.format(type(self.file_path), self.file_path))
