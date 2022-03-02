import os
import pytest

from intervaltree import IntervalTree
from cte import amplicon_scheme

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_scheme")


def test_load_primers_file():
    infile = os.path.join(data_dir, "load_primers_file.tsv")
    got_primers, got_amplicons = amplicon_scheme.AmpliconScheme._load_primers_file(
        infile
    )
    expect_primers = IntervalTree()
    expect_primers[30:36] = "amp1_left"
    expect_primers[29:31] = "amp1_left_alt"
    expect_primers[100:104] = "amp1_right"
    expect_primers[98:104] = "amp2_left"
    expect_primers[150:157] = "amp2_right"
    expect_primers[145:151] = "amp3_left"
    expect_primers[190:195] = "amp3_right"
    expect_amplicons = {
        "amp1": {
            "previous": None,
            "next": "amp2",
            "start": 29,
            "end": 103,
            "left_primers": [(30, 35), (29, 30)],
            "right_primers": [(100, 103)],
        },
        "amp2": {
            "previous": "amp1",
            "next": "amp3",
            "start": 98,
            "end": 156,
            "left_primers": [(98, 103)],
            "right_primers": [(150, 156)],
        },
        "amp3": {
            "previous": "amp2",
            "next": None,
            "start": 145,
            "end": 194,
            "left_primers": [(145, 150)],
            "right_primers": [(190, 194)],
        },
    }
    assert got_primers == expect_primers
    assert got_amplicons == expect_amplicons


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
    assert scheme.amplicon_tree == expect_amplicons
    assert not scheme.in_primer(27)
    assert not scheme.in_primer(28)
    assert scheme.in_primer(29)
    assert scheme.in_primer(32)
    assert not scheme.in_primer(33)
    assert not scheme.overlaps_primer(28, 28)
    assert scheme.overlaps_primer(28, 29)
    assert scheme.overlaps_primer(28, 32)
    assert not scheme.overlaps_primer(33, 34)
    assert scheme.overlaps_primer(33, 130)
    assert scheme.overlapping_amplicons(1, 2) is None
    assert scheme.overlapping_amplicons(1, 28) is None
    assert scheme.overlapping_amplicons(1, 29) == ["amp1"]
    assert scheme.overlapping_amplicons(1, 130) == ["amp1", "amp2"]
    assert scheme.overlapping_amplicons(145, 147) is None

    assert scheme.amp_coords("amp1") == (29, 43)
    assert scheme.amp_coords("amp2") == (130, 144)
    assert scheme.inner_amp_coords("amp1") == (33, 129)
    assert scheme.inner_amp_coords("amp2") == (44, 141)
