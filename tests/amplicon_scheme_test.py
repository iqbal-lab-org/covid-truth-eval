import os
import pytest

from intervaltree import IntervalTree, Interval
from cte import amplicon_scheme

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_scheme")


def test_amplicon_scheme():
    infile = os.path.join(data_dir, "load_viridian_workflow_primers_tsv.tsv")
    scheme = amplicon_scheme.AmpliconScheme(infile)
    expect_primers = IntervalTree()
    expect_primers[30:33] = "amp1_left"
    expect_primers[29:31] = "amp1_left_alt"
    expect_primers[42:44] = "amp1_right"
    expect_primers[130:133] = "amp2_left"
    expect_primers[142:144] = "amp2_right"
    expect_primers[143:145] = "amp2_right_alt"

    expect_amplicons = IntervalTree()
    expect_amplicons[29:44] = "amp1"
    expect_amplicons[130:145] = "amp2"

    assert scheme.primers == expect_primers
    assert scheme.first_primer_start == 29
    assert scheme.last_primer_end == 144
    assert scheme.amplicons == expect_amplicons
    assert not scheme.overlaps_primer(28, 28)
    assert scheme.overlaps_primer(28, 29)
    assert scheme.overlaps_primer(28, 32)
    assert not scheme.overlaps_primer(33, 34)
    assert scheme.overlaps_primer(33, 130)
    assert scheme.overlapping_amplicons(1, 2) is None
    assert scheme.overlapping_amplicons(1, 28) is None
    amp1 = {"name": "amp1", "start": 29, "end": 43}
    amp2 = {"name": "amp2", "start": 130, "end": 144}
    assert scheme.overlapping_amplicons(1, 29) == [amp1]
    assert scheme.overlapping_amplicons(1, 130) == [amp1, amp2]
    assert scheme.overlapping_amplicons(145, 147) is None
