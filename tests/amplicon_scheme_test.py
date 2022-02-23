import os
import pytest

from cte import primers

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "primers")


def test_load_viridian_workflow_primers_tsv():
    infile = os.path.join(data_dir, "load_viridian_workflow_primers_tsv.tsv")
    got = primers.load_viridian_workflow_primers_tsv(infile)
    expect_lookup = {
        29: {("amp1", "amp1_left_alt")},
        30: {("amp1", "amp1_left"), ("amp1", "amp1_left_alt")},
        31: {("amp1", "amp1_left")},
        32: {("amp1", "amp1_left")},
        42: {("amp1", "amp1_right")},
        43: {("amp1", "amp1_right")},
        130: {("amp2", "amp2_left")},
        131: {("amp2", "amp2_left")},
        132: {("amp2", "amp2_left")},
        142: {("amp2", "amp2_right")},
        143: {("amp2", "amp2_right"), ("amp2", "amp2_right_alt")},
        144: {("amp2", "amp2_right_alt")},
    }
    expect_start = 29
    expect_end = 144
    assert got == (expect_lookup, expect_start, expect_end)
