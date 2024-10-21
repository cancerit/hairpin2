# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from hairpin2 import main as hp2
import pysam
import pytest
import factory
import factory.random
from faker import Faker
from faker_biology.bioseq import Bioseq
import random

factory.random.reseed_random(2501)
random.seed(2501)


# smoke test validate_read
# ----
class ExtendedBioProvider(Bioseq):
    def quality_string(self, length):
        if length < 1:
            raise ValueError('length must be geater than 1')
        allowed_chars = [chr(i) for i in range(33, 75)]
        return ''.join([random.choice(allowed_chars) for _ in range(length)])

    def cigar_string(self, length):
        if length < 1:
            raise ValueError('length must be greater than 1')
        opchars = 'MIDNSHP=XB'
        opchars_noclip = 'MIDNP=XB'
        bound = 200 if length > 200 else length
        cig_op_lengths = []
        while (bound > 0):
            oplen = random.randint(1, bound)
            cig_op_lengths.append(random.randint(1, oplen))
            cig_sum = sum(cig_op_lengths)
            bound = 200 if length - cig_sum > 200 else length - cig_sum
        cig_op_lengths[-1] = cig_op_lengths[-1] - (cig_sum - length)
        cig_ops = []
        last_opchar = ''
        # first and last op can be S or H, but not others
        # first op H last op S, i.e. only clipping ops, seg faults pysam
        # reads with only clipping ops seem to segfault pysam... report bug
        if len(cig_op_lengths) == 1:
            cig_ops.append(random.choice(opchars_noclip))
        else:
            cig_ops.append(random.choice(opchars))
        for _ in range(max([len(cig_op_lengths) - 2, 0])):
            iter_opchars = opchars_noclip.replace(last_opchar, '')
            cig_ops.append(random.choice(iter_opchars))
            last_opchar = cig_ops[-1]
        if len(cig_ops) != 1:
            cig_ops.append(random.choice(opchars_noclip if cig_ops[-1] in ['H', 'S'] else opchars))
        return ''.join([str(x) for pair in zip(cig_op_lengths, cig_ops) for x in pair])


fake = Faker()
fake.add_provider(ExtendedBioProvider)


class AlignedSegmentWrapper:
    def __init__(self, query_name, query_sequence, query_qualities, flag, reference_id, reference_start, next_reference_start, mapping_quality, cigarstring, mc):
        self.segment = pysam.AlignedSegment()
        self.segment.query_name = query_name
        self.segment.query_sequence = query_sequence
        self.segment.query_qualities = pysam.qualitystring_to_array(query_qualities)
        self.segment.flag = flag
        self.segment.reference_id = reference_id
        self.segment.reference_start = reference_start
        self.segment.next_reference_start = next_reference_start
        self.segment.mapping_quality = mapping_quality
        self.segment.cigarstring = cigarstring
        self.segment.set_tag('MC', mc)


class ReadFactory(factory.Factory):
    class Meta:
        model = AlignedSegmentWrapper

    query_name = 'read1'  # should one assume pysam handles all bizarre query names gracefully? I am...
    query_sequence = factory.LazyAttribute(lambda _: fake.dna(length=random.randint(50, 200)))
    query_qualities = factory.LazyAttribute(lambda o: fake.quality_string(length=len(o.query_sequence)))
    flag = factory.LazyAttribute(lambda _: random.getrandbits(16))
    reference_id = 0
    reference_start = factory.LazyAttribute(lambda _: random.randint(1, 300000000))
    next_reference_start = factory.LazyAttribute(lambda o: o.reference_start - random.randint(-700, 700))
    mapping_quality = factory.LazyAttribute(lambda _: random.randint(0, 255))
    cigarstring = factory.LazyAttribute(lambda o: fake.cigar_string(length=len(o.query_sequence)))
    mc = factory.LazyAttribute(lambda _: fake.cigar_string(length=random.randint(50, 200)))


@pytest.mark.dev
@pytest.mark.parametrize("test_read", [ReadFactory().segment for _ in range(1000)])
def test_smoke(test_read):
    mut_pos = random.randint(1, len(test_read.query_sequence) - 1)
    start = test_read.reference_start + mut_pos
    alt = random.choices([test_read.query_sequence[mut_pos:mut_pos + random.randint(1, 3)], '.'], cum_weights=[66, 100])[0]
    if alt == '.':
        mut_type_str = 'D'
        stop = start + random.randint(1, 3)
    elif len(alt) == 1:
        mut_type_str = random.choice(['S', 'I'])
        stop = start + 1
    else:
        mut_type_str = random.choice(['S', 'I'])
        stop = start + 1 if mut_type_str == 'I' else start + len(alt)
    vflag = hp2.validate_read(test_read,
                              vcf_start=start,
                              vcf_stop=stop,
                              alt=alt,
                              mut_type=mut_type_str,
                              min_mapqual=11,
                              min_clipqual=35,
                              min_basequal=25)
    print(format(vflag, '010b'))
